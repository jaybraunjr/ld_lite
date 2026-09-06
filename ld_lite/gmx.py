"""The two GROMACS input files a built LD needs: ``topol.top`` and ``index.ndx``."""

import warnings

import MDAnalysis as mda

#: include order for the .top; only residues actually present are written
ITP_ORDER = ["POPC", "DOPE", "SAPI", "POPE", "POPG", "POPS", "PSM", "CHL1",
             "CHYO", "TRIO", "SOD", "CLA", "POT", "TIP3"]

#: itp filenames that don't follow the "<RESNAME>.itp" convention
ITP_FILENAME = {"CHYO": "chyo.itp"}

#: index.ndx buckets, mapped onto the physical layers of the droplet.
#: CHYO is a *membrane* lipid here: ld_lite only ever inserts TRIO, so any
#: cholesteryl ester in the system came from the bilayer's leaflets.
CORE = ["TRIO"]
MEMBRANE = ["POPC", "DOPE", "SAPI", "POPE", "POPG", "POPS", "PSM", "CHL1", "CHYO"]
SOLVENT = ["TIP3", "SOL", "SOD", "CLA", "POT", "K", "NA"]


def molecule_counts(universe):
    """Ordered ``(resname, count)`` runs, in file order.

    GROMACS reads ``[ molecules ]`` positionally, so consecutive runs of the
    same residue are collapsed rather than totalled per name.
    """
    runs = []
    for res in universe.residues:
        if runs and runs[-1][0] == res.resname:
            runs[-1][1] += 1
        else:
            runs.append([res.resname, 1])
    return [(rn, n) for rn, n in runs]


def write_top(universe, path, includes=None,
              system_name="Lipid droplet (bilayer + TRIO core)", known=()):
    """Write a ``.top`` with force-field includes and the ``[ molecules ]`` block.

    ``known`` names extra residues to accept beyond :data:`ITP_ORDER`; each is
    included as ``toppar/<RESNAME>.itp`` after the built-in ones, so a lipid the
    package has never heard of needs no edit to this file.
    """
    if includes is None:
        resnames = set(universe.residues.resnames)
        order = list(ITP_ORDER) + [r for r in known if r not in ITP_ORDER]
        present = [r for r in order if r in resnames]
        unknown = resnames - set(order)
        if unknown:
            warnings.warn(
                f"no #include guessed for {sorted(unknown)} -- add them to the "
                f"topology by hand, or pass includes=[...]", RuntimeWarning)
        includes = ["toppar/charmm36-mar2019.ff/forcefield.itp"] + \
                   [f"toppar/{ITP_FILENAME.get(r, r + '.itp')}" for r in present]
    with open(path, "w") as fh:
        for inc in includes:
            fh.write(f'#include "{inc}"\n')
        fh.write(f"\n[ system ]\n{system_name}\n\n")
        fh.write("[ molecules ]\n; Compound        #mols\n")
        for rn, n in molecule_counts(universe):
            fh.write(f"{rn:<16} {n}\n")
    return path


def write_index(structure, path, extra_core=(), extra_membrane=(), extra_solvent=()):
    """Write the ``index.ndx`` groups the CHARMM-GUI protocol expects.

    For a droplet the thermostat groups map onto the layers: ``SOLU`` is the
    TRIO core, ``MEMB`` the phospholipids plus any surface CHYO, ``SOLV`` the
    water and ions.  That keeps all three non-empty without a protein in the
    system.  Raises if any residue would be left out of all three.

    The ``extra_*`` arguments add residues to a group for this call only, so an
    unrecognised lipid can be placed from the command line rather than by
    editing the lists above.
    """
    u = structure if isinstance(structure, mda.Universe) else mda.Universe(structure)
    present = set(u.residues.resnames)
    core = [r for r in list(CORE) + list(extra_core) if r in present]
    memb = [r for r in list(MEMBRANE) + list(extra_membrane) if r in present]
    solv = [r for r in list(SOLVENT) + list(extra_solvent) if r in present]

    # every atom must land in exactly one thermostat group, or grompp rejects
    # the run with "atom ... does not belong to any T-coupling group"
    unassigned = present - set(core) - set(memb) - set(solv)
    if unassigned:
        raise ValueError(
            f"residues {sorted(unassigned)} fall into no index group -- declare "
            f"them with --membrane/--solvent/--core (or the extra_* arguments), "
            f"otherwise grompp will reject the system for having atoms outside "
            f"tc_grps")

    def sel(names):
        return "resname " + " ".join(names) if names else "resid -1"

    groups = {
        "SOLU": sel(core),
        "MEMB": sel(memb),
        "SOLV": sel(solv),
        "SOLU_MEMB": sel(core + memb),
        "SYSTEM": "all",
    }
    with open(path, "w") as fh:
        for name, selection in groups.items():
            idx = u.select_atoms(selection).indices + 1        # GROMACS is 1-based
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(map(str, idx[i:i + 15])) + "\n")
    return {name: int(len(u.select_atoms(s))) for name, s in groups.items()}
