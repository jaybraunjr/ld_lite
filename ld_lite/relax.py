"""Rigid-body clash relaxation for the freshly placed core.

Molecules dropped onto a grid leave hard atomic overlaps that would blow up an
MD minimiser.  This nudges whole molecules apart -- translating each as a rigid
body, never distorting its conformation -- until the closest contact clears
``d_min``.

The bilayer is already equilibrated, so it is passed in as ``obstacles``: the
core is pushed away from the leaflets and the water, but nothing in the bilayer
moves.  Relaxing the lipids too would shove them sideways into water molecules
that stay put, which creates more clashes than it removes.
"""

import numpy as np
from MDAnalysis.lib.distances import capped_distance


def _push(a, b, dists, d_min, step, box):
    """Displacement that separates each contacting pair, damped by ``step``."""
    overlap = d_min - np.maximum(dists, 1e-3)
    direction = a - b
    direction -= box[:3] * np.round(direction / box[:3])          # minimum image
    direction /= np.linalg.norm(direction, axis=1, keepdims=True) + 1e-9
    return (overlap * step)[:, None] * direction


def relax_clashes(molecules, box, obstacles=None, d_min=1.6, iterations=60,
                  step=0.4, max_move=1.5, verbose=False):
    """Separate clashing molecules in place, and report the closest contact.

    Parameters
    ----------
    molecules : list[(resname, names, coords)]
        The molecules free to move; ``coords`` arrays are modified in place.
    box : array-like (6,)
        MDAnalysis box (a, b, c, alpha, beta, gamma) for minimum-image contacts.
    obstacles : ndarray (M, 3), optional
        Fixed atoms that repel the molecules but never move themselves.
    d_min : float
        Target minimum inter-molecular atom-atom distance, angstrom.
    step : float
        Fraction of each overlap moved per sweep.
    max_move : float
        Cap on how far one molecule may travel in a single sweep, angstrom.

    Returns
    -------
    float
        The closest inter-molecular contact remaining, in angstrom.
    """
    n = len(molecules)
    sizes = np.array([len(m[1]) for m in molecules])
    mol_of_atom = np.repeat(np.arange(n), sizes)
    box = np.asarray(box, dtype=np.float32)
    if obstacles is not None:
        obstacles = np.ascontiguousarray(obstacles, dtype=np.float32)

    closest = float("inf")
    for it in range(iterations):
        coords = np.vstack([m[2] for m in molecules]).astype(np.float32)
        shift = np.zeros((n, 3))
        n_contacts = 0
        closest = float("inf")

        # --- molecule against molecule ---
        pairs, dists = capped_distance(coords, coords, max_cutoff=d_min, box=box,
                                       return_distances=True)
        # genuine inter-molecular contacts, counted once
        mask = mol_of_atom[pairs[:, 0]] != mol_of_atom[pairs[:, 1]]
        mask &= pairs[:, 0] < pairs[:, 1]
        pairs, dists = pairs[mask], dists[mask]
        if len(pairs):
            push = _push(coords[pairs[:, 0]], coords[pairs[:, 1]], dists, d_min, step, box)
            # sum rather than average, so a deeply buried molecule can still escape
            np.add.at(shift, mol_of_atom[pairs[:, 0]], push)
            np.add.at(shift, mol_of_atom[pairs[:, 1]], -push)
            n_contacts += len(pairs)
            closest = min(closest, float(dists.min()))

        # --- molecule against the fixed bilayer, pushing only our side ---
        if obstacles is not None and len(obstacles):
            opairs, odists = capped_distance(coords, obstacles, max_cutoff=d_min,
                                             box=box, return_distances=True)
            if len(opairs):
                push = _push(coords[opairs[:, 0]], obstacles[opairs[:, 1]],
                             odists, d_min, step, box)
                np.add.at(shift, mol_of_atom[opairs[:, 0]], push)
                n_contacts += len(opairs)
                closest = min(closest, float(odists.min()))

        if n_contacts == 0:
            if verbose:
                print(f"  relaxed: clash-free after {it} sweeps")
            return d_min

        mag = np.linalg.norm(shift, axis=1)
        over = mag > max_move
        shift[over] *= (max_move / mag[over])[:, None]

        for mi in np.nonzero(mag > 1e-6)[0]:
            molecules[mi][2][:] = molecules[mi][2] + shift[mi]

        if verbose and it % 10 == 0:
            print(f"  sweep {it}: {n_contacts} contacts under {d_min} A, "
                  f"closest {closest:.2f} A")
    return closest
