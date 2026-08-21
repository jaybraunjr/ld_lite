"""Turn a phospholipid bilayer into a lipid droplet by inserting a TRIO core.

Start from an equilibrated, solvated, ionized CHARMM-GUI bilayer, pry the two
leaflets apart along z, and fill the vacated midplane with triolein.  The
bilayer's own water and ions ride along untouched -- they sit outside the
membrane -- so only the hydrophobic interior changes.  The result is the LD
"trilayer": PL monolayer / TRIO core / PL monolayer.

    bilayer                          lipid droplet
  ~~~~~~~~~~~~  heads              ~~~~~~~~~~~~   monolayer
  ||||||||||||  tails      -->     ||||||||||||
  ||||||||||||  tails              ############   <- inserted TRIO core
  ~~~~~~~~~~~~  heads              ||||||||||||
                                   ~~~~~~~~~~~~   monolayer
"""

import os

import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree

from .templates import extract_template
from .relax import relax_clashes

#: packaged single-molecule TRIO template, so no external structure is needed
DEFAULT_CORE_SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "data", "core_templates.gro")

#: approximate molecular volume of triolein (885 g/mol at ~0.91 g/cm^3)
TRIO_VOLUME_A3 = 1630.0


def _midplane_z(universe, head_name="P"):
    """Membrane mid-plane z, taken as the mean height of all head-group atoms."""
    heads = universe.select_atoms(f"name {head_name}")
    if not len(heads):
        raise ValueError(
            f"no '{head_name}' atoms found -- is this a phospholipid bilayer?")
    return float(heads.positions[:, 2].mean())


def make_whole(universe):
    """Reassemble residues that are split across the periodic boundary.

    A ``.gro`` carries no bonds, so :meth:`AtomGroup.unwrap` is unavailable and
    ``wrap(compound="residues")`` is no help -- it translates each residue by
    its centre of mass, which for an already-split residue is a meaningless
    point midway across the box.  Instead every atom is moved to its minimum
    image relative to the first atom of its own residue, which reassembles any
    molecule smaller than half the box without needing connectivity.

    This matters twice over: a split residue has a bogus centre of mass, which
    would put it in the wrong leaflet, and it makes ``grompp`` abort with
    "largest distance between excluded atoms ... larger than the cut-off".

    Returns the number of residues that had to be repaired.
    """
    box = universe.dimensions[:3].astype(float)
    if np.any(box <= 0):
        raise ValueError("the bilayer has no box -- cannot resolve periodicity")

    pos = universe.atoms.positions.astype(float)
    resindices = universe.atoms.resindices
    # first atom of each residue, as a per-atom reference point
    _, first = np.unique(resindices, return_index=True)
    ref = pos[first][resindices]

    shift = box * np.round((pos - ref) / box)
    moved = np.abs(shift).sum(axis=1) > 0
    if moved.any():
        universe.atoms.positions = pos - shift
    return int(np.unique(resindices[moved]).size)


def build_ld(bilayer, n_trio, core_source=None, core_thickness=None,
             packing_fraction=0.25, head_name="P", seed=0, relax=True,
             d_min=1.6):
    """Insert a TRIO core between a bilayer's leaflets and return the LD.

    Parameters
    ----------
    bilayer : str or mda.Universe
        A solvated phospholipid bilayer, e.g. a CHARMM-GUI ``.gro``.
    n_trio : int
        Number of triolein molecules to place in the core.
    core_source : str or mda.Universe, optional
        Structure to take the TRIO conformation from.  Defaults to the
        packaged template.
    core_thickness : float, optional
        Thickness of the inserted slab in angstrom.  Derived from ``n_trio``
        and the box area via ``packing_fraction`` when left unset.
    packing_fraction : float
        Target initial fill of the core, 0-1.  This is the single biggest lever
        on how clean the structure comes out: a looser core leaves fewer
        overlaps but a taller box for NPT to compress back down.
    head_name : str
        Atom name marking phospholipid head groups (``P`` for PC/PE/PI/PS/PG).
    relax : bool
        Run rigid-body clash relaxation on the core and lipids afterwards.

    Returns
    -------
    mda.Universe
        The LD system, with the box grown along z to fit the core.
    """
    n_trio = int(n_trio)
    if n_trio < 1:
        raise ValueError("n_trio must be at least 1")

    u = bilayer if isinstance(bilayer, mda.Universe) else mda.Universe(bilayer)
    # reassemble anything split across the boundary, *then* put whole residues
    # back in the box -- in that order, so the leaflet split sees real centres
    repaired = make_whole(u)
    if repaired:
        print(f"[ld_lite] repaired {repaired} residues split across the periodic boundary")
    u.atoms.wrap(compound="residues")
    box = u.dimensions.copy()
    centre = _midplane_z(u, head_name)

    if core_thickness is None:
        core_thickness = (n_trio * TRIO_VOLUME_A3) / (box[0] * box[1] * packing_fraction)
    t = float(core_thickness)

    # --- open the gap: every residue moves with whichever half it belongs to ---
    res_com_z = np.array([r.atoms.center_of_mass()[2] for r in u.residues])
    top = res_com_z >= centre
    u.residues[top].atoms.positions += np.array([0.0, 0.0, +t / 2.0])
    u.residues[~top].atoms.positions += np.array([0.0, 0.0, -t / 2.0])

    tpl = extract_template(core_source or DEFAULT_CORE_SOURCE, "TRIO")
    core_mols = _place_core(u, tpl, n_trio, box, centre, t, np.random.default_rng(seed))

    ld = _merge(u, core_mols, box, t)
    # recentre so the core sits at the new box midplane, water symmetric about it
    ld.atoms.translate([0.0, 0.0, ld.dimensions[2] / 2.0 - centre])

    if relax:
        # only the core moves; the bilayer is already equilibrated, so it acts
        # as a fixed obstacle that the core is pushed away from
        core = ld.select_atoms("resname TRIO")
        mols = [("TRIO", list(r.atoms.names), r.atoms.positions.copy())
                for r in core.residues]
        closest = relax_clashes(mols, box=ld.dimensions, d_min=d_min,
                                obstacles=ld.select_atoms("not resname TRIO").positions)
        ld.atoms[core.indices].positions = np.vstack([m[2] for m in mols])
        print(f"[ld_lite] closest contact after relaxation: {closest:.2f} A")
    return ld


def _place_core(u, tpl, n_trio, box, centre, t, rng, tries=100, accept=2.0):
    """Drop ``n_trio`` randomly oriented TRIO onto a jittered 3-D grid in the gap.

    Each molecule is tried in several orientations and the one that sits
    furthest from everything already placed wins.  Rigid relaxation cannot
    un-thread interlocked molecules, so the trick is to never interlock them:
    the occupancy is seeded with the leaflet atoms bordering the gap, and the
    core avoids the tails as well as itself.
    """
    span = tpl.span
    cell = 0.85 * float(max(span[0], span[1]))
    nx = max(1, int(box[0] // cell))
    ny = max(1, int(box[1] // cell))
    nz = max(1, int(np.ceil(n_trio / (nx * ny))))
    sx, sy = box[0] / nx, box[1] / ny

    # inset the centres so molecules don't poke far past the core interfaces
    inset = min(0.30 * t, 0.25 * float(span[2]))
    z_lo, z_hi = centre - t / 2 + inset, centre + t / 2 - inset
    sz = (z_hi - z_lo) / nz
    jit = 0.2 * min(sx, sy, max(sz, 1.0))

    slots = [(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
    rng.shuffle(slots)

    border = (u.atoms.positions[:, 2] > z_lo - 4) & (u.atoms.positions[:, 2] < z_hi + 4)
    placed = u.atoms.positions[border].astype(np.float32)

    core_mols = []
    for i, j, k in slots[:n_trio]:
        base = np.array([(i + 0.5) * sx, (j + 0.5) * sy, z_lo + (k + 0.5) * sz])
        tree = cKDTree(placed) if len(placed) else None
        best, best_min = None, -1.0
        for _ in range(tries):
            coords = tpl.randomly_rotated(rng) + base + rng.uniform(-jit, jit, 3)
            if tree is None:
                best = coords
                break
            mind = float(tree.query(coords, k=1)[0].min())
            if mind > best_min:
                best_min, best = mind, coords
            if mind >= accept:
                break
        core_mols.append(("TRIO", list(tpl.names), best))
        placed = np.vstack([placed, best]).astype(np.float32)
    return core_mols


def _merge(u, core_mols, box, t):
    """Combine the pried-apart bilayer with the inserted core molecules."""
    core_counts = [len(names) for _, names, _ in core_mols]
    n_res = u.residues.n_residues + len(core_mols)
    counts = [r.atoms.n_atoms for r in u.residues] + core_counts

    new = mda.Universe.empty(
        n_atoms=u.atoms.n_atoms + sum(core_counts),
        n_residues=n_res,
        atom_resindex=np.repeat(np.arange(n_res), counts),
        trajectory=True,
    )
    new.add_TopologyAttr(
        "names", list(u.atoms.names) + [nm for _, names, _ in core_mols for nm in names])
    new.add_TopologyAttr(
        "resnames", list(u.residues.resnames) + [rn for rn, _, _ in core_mols])
    new.add_TopologyAttr("resids", np.arange(1, n_res + 1))
    new.atoms.positions = np.vstack([u.atoms.positions] + [c for _, _, c in core_mols])

    box = box.copy()
    box[2] += t                                  # grow z to fit the core
    new.dimensions = box
    return new
