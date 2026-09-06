"""Fast tests on a synthetic bilayer -- no large files and no GROMACS needed."""

import re

import numpy as np
import MDAnalysis as mda
import pytest

from ld_lite import (build_ld, extract_template, make_whole, molecule_counts,
                     write_top, write_index)
from ld_lite import equil
from ld_lite.build import DEFAULT_CORE_SOURCE

BOX = (60.0, 60.0, 80.0)
MIDPLANE = 40.0


def synthetic_bilayer(n_side=4, n_water=20, n_chyo=0):
    """A 4-atom-per-lipid bilayer with water, enough for the geometry to work on.

    ``n_chyo`` adds cholesteryl-ester residues to the leaflets; they carry no
    ``P`` atom, exactly like the real thing.
    """
    lipid_names = ["P", "C1", "C2", "C3"]
    lipid_masses = [31.0, 12.0, 12.0, 12.0]
    water_names = ["OH2", "H1", "H2"]
    water_masses = [16.0, 1.0, 1.0]

    resnames, names, masses, coords = [], [], [], []
    pitch = BOX[0] / n_side
    for leaflet, (head_z, direction) in enumerate([(55.0, -1), (25.0, +1)]):
        for i in range(n_side):
            for j in range(n_side):
                x, y = (i + 0.5) * pitch, (j + 0.5) * pitch
                resnames.append("POPC")
                names.extend(lipid_names)
                masses.extend(lipid_masses)
                for k in range(4):
                    coords.append([x, y, head_z + direction * 5.0 * k])

    rng = np.random.default_rng(0)
    for w in range(n_water):
        z = 10.0 if w % 2 else 70.0
        x, y = rng.uniform(0, BOX[0], 2)
        resnames.append("TIP3")
        names.extend(water_names)
        masses.extend(water_masses)
        coords.extend([[x, y, z], [x + 0.9, y, z], [x, y + 0.9, z]])

    for c in range(n_chyo):
        # sits in a leaflet, head-group free
        z = 50.0 if c % 2 else 30.0
        x, y = rng.uniform(0, BOX[0], 2)
        resnames.append("CHYO")
        names.extend(["C1", "C2", "C3", "C4"])
        masses.extend([12.0] * 4)
        coords.extend([[x, y, z + k] for k in range(4)])

    counts = [4] * (2 * n_side * n_side) + [3] * n_water + [4] * n_chyo
    u = mda.Universe.empty(
        n_atoms=sum(counts), n_residues=len(counts),
        atom_resindex=np.repeat(np.arange(len(counts)), counts), trajectory=True)
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", np.arange(1, len(counts) + 1))
    u.add_TopologyAttr("masses", masses)
    u.atoms.positions = np.array(coords, dtype=float)
    u.dimensions = [*BOX, 90.0, 90.0, 90.0]
    return u


@pytest.fixture(scope="module")
def ld():
    return build_ld(synthetic_bilayer(), n_trio=4, core_thickness=30.0, seed=1)


def test_packaged_trio_template_is_centred():
    tpl = extract_template(DEFAULT_CORE_SOURCE, "TRIO")
    assert tpl.n_atoms > 100
    assert np.allclose(tpl.positions.mean(axis=0), 0.0, atol=1e-3)


def test_missing_residue_is_reported():
    with pytest.raises(ValueError, match="NOPE"):
        extract_template(DEFAULT_CORE_SOURCE, "NOPE")


def test_bilayer_without_head_atoms_is_rejected():
    u = synthetic_bilayer()
    with pytest.raises(ValueError, match="bilayer"):
        build_ld(u, n_trio=2, head_name="ZZ")


def test_core_is_inserted_and_box_grows(ld):
    assert len(ld.select_atoms("resname TRIO").residues) == 4
    assert len(ld.select_atoms("resname POPC").residues) == 32
    assert ld.dimensions[2] == pytest.approx(BOX[2] + 30.0)
    assert ld.dimensions[0] == pytest.approx(BOX[0])


def test_core_sits_between_the_leaflets(ld):
    core_z = ld.select_atoms("resname TRIO").positions[:, 2]
    heads_z = ld.select_atoms("name P").positions[:, 2]
    # the core straddles the midplane, the heads stay outside it
    assert heads_z.min() < core_z.mean() < heads_z.max()
    assert core_z.mean() == pytest.approx(ld.dimensions[2] / 2.0, abs=8.0)


def _split_one_residue(u, resid_index=0):
    """Push half of one residue across the x boundary, as trjconv output does."""
    res = u.residues[resid_index]
    pos = u.atoms.positions.copy()
    idx = res.atoms.ix[len(res.atoms) // 2:]
    pos[idx, 0] += BOX[0]
    u.atoms.positions = pos
    return idx


def test_make_whole_reassembles_a_split_residue():
    u = synthetic_bilayer()
    idx = _split_one_residue(u)
    span_before = u.residues[0].atoms.positions[:, 0].ptp()
    assert span_before > BOX[0] / 2          # genuinely split

    repaired = make_whole(u)
    assert repaired == 1
    assert u.residues[0].atoms.positions[:, 0].ptp() < BOX[0] / 2
    # every other residue is left alone
    assert make_whole(u) == 0


def test_split_residue_does_not_survive_into_the_build():
    """grompp aborts if any molecule spans more than the cut-off, so the
    builder must reassemble the input rather than propagate the break."""
    u = synthetic_bilayer()
    _split_one_residue(u)
    ld = build_ld(u, n_trio=4, core_thickness=30.0, seed=1)

    # No bilayer residue may span more than half the box along any axis; a
    # longer one is still split. (Only the input is checked -- a single TRIO is
    # ~33 A end to end, which legitimately exceeds half of this toy 60 A box.)
    half = ld.dimensions[:3] / 2.0
    bilayer_res = [r for r in ld.residues if r.resname != "TRIO"]
    worst = max(float((r.atoms.positions.ptp(axis=0) / half).max()) for r in bilayer_res)
    assert worst < 1.0


def test_molecule_counts_collapse_runs():
    u = synthetic_bilayer(n_side=2, n_water=3)
    assert molecule_counts(u) == [("POPC", 8), ("TIP3", 3)]


def test_top_and_index_are_written(ld, tmp_path):
    top = write_top(ld, tmp_path / "topol.top")
    text = open(top).read()
    assert re.search(r"^TRIO\s+4$", text, re.M)
    assert re.search(r"^POPC\s+32$", text, re.M)
    assert '#include "toppar/POPC.itp"' in text

    counts = write_index(ld, tmp_path / "index.ndx")
    assert counts["SOLU"] == len(ld.select_atoms("resname TRIO"))
    assert counts["MEMB"] == len(ld.select_atoms("resname POPC"))
    assert counts["SOLV"] == len(ld.select_atoms("resname TIP3"))
    assert counts["SOLU_MEMB"] == counts["SOLU"] + counts["MEMB"]
    assert counts["SYSTEM"] == ld.atoms.n_atoms
    assert "[ SOLU_MEMB ]" in open(tmp_path / "index.ndx").read()


def test_chyo_survives_the_build_and_lands_in_memb(tmp_path):
    """CHYO in the leaflets is a membrane lipid, not part of the inserted core."""
    ld = build_ld(synthetic_bilayer(n_chyo=6), n_trio=4, core_thickness=30.0, seed=1)
    assert len(ld.select_atoms("resname CHYO").residues) == 6

    counts = write_index(ld, tmp_path / "index.ndx")
    chyo = len(ld.select_atoms("resname CHYO"))
    assert counts["SOLU"] == len(ld.select_atoms("resname TRIO"))   # core is TRIO only
    assert counts["MEMB"] == len(ld.select_atoms("resname POPC")) + chyo
    # every atom must be covered by the three thermostat groups
    assert counts["SOLU"] + counts["MEMB"] + counts["SOLV"] == counts["SYSTEM"]


def test_chyo_include_uses_the_lowercase_filename(tmp_path):
    ld = build_ld(synthetic_bilayer(n_chyo=6), n_trio=4, core_thickness=30.0, seed=1)
    text = write_top(ld, tmp_path / "topol.top") and (tmp_path / "topol.top").read_text()
    assert '#include "toppar/chyo.itp"' in text
    assert '#include "toppar/CHYO.itp"' not in text
    assert re.search(r"^CHYO\s+6$", text, re.M)


def test_residue_in_no_index_group_is_rejected(tmp_path):
    u = synthetic_bilayer(n_side=2, n_water=2)
    u.residues[0].resname = "XXXX"
    with pytest.raises(ValueError, match="XXXX"):
        write_index(u, tmp_path / "index.ndx")


def test_unknown_lipid_can_be_declared_as_membrane(tmp_path):
    """A lipid the package has never heard of must be usable without a source edit."""
    u = synthetic_bilayer(n_side=2, n_water=2)
    u.residues[0].resname = "DPPC"

    counts = write_index(u, tmp_path / "index.ndx", extra_membrane=["DPPC"])
    assert counts["MEMB"] == len(u.select_atoms("resname POPC DPPC"))
    assert counts["SOLU"] + counts["MEMB"] + counts["SOLV"] == counts["SYSTEM"]

    write_top(u, tmp_path / "topol.top", known=["DPPC"])
    text = (tmp_path / "topol.top").read_text()
    assert '#include "toppar/DPPC.itp"' in text
    assert re.search(r"^DPPC\s+1$", text, re.M)


def test_declared_residue_can_go_to_solvent_or_core(tmp_path):
    u = synthetic_bilayer(n_side=2, n_water=3)
    u.residues[-1].resname = "URE"                      # a cosolvent, say

    counts = write_index(u, tmp_path / "index.ndx", extra_solvent=["URE"])
    assert counts["SOLV"] == len(u.select_atoms("resname TIP3 URE"))
    # the same residue routed to the core instead lands in SOLU, not SOLV
    counts = write_index(u, tmp_path / "index.ndx", extra_core=["URE"])
    assert counts["SOLU"] == len(u.select_atoms("resname URE"))


def test_undeclared_lipid_error_names_the_flags(tmp_path):
    u = synthetic_bilayer(n_side=2, n_water=2)
    u.residues[0].resname = "DPPC"
    with pytest.raises(ValueError, match="--membrane"):
        write_index(u, tmp_path / "index.ndx")


def test_mdps_are_packaged(tmp_path):
    equil.stage_mdps(tmp_path)
    staged = sorted(p.name for p in tmp_path.iterdir())
    assert staged == sorted([equil.MINIMIZATION] + equil.EQUILIBRATION + [equil.PRODUCTION])
    assert "integrator              = steep" in (tmp_path / equil.MINIMIZATION).read_text()


def test_equil_refuses_an_unbuilt_directory(tmp_path):
    with pytest.raises(FileNotFoundError, match="ld_lite build"):
        equil.run(str(tmp_path))


def test_equil_requires_toppar(ld, tmp_path):
    for name in ("system.gro", "ref.gro"):
        ld.atoms.write(tmp_path / name)
    write_top(ld, tmp_path / "topol.top")
    write_index(ld, tmp_path / "index.ndx")
    with pytest.raises(FileNotFoundError, match="toppar"):
        equil.run(str(tmp_path))
