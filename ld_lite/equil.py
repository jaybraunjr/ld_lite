"""Run the CHARMM-GUI minimization + equilibration protocol on a staged run.

``step6.0`` is energy minimization; ``step6.1``-``step6.6`` walk the position
restraints down from stiff to free while the timestep grows from 1 to 2 fs and
the barostat switches on.  ``step7`` is unrestrained NPT production, run only
when asked for.  Every stage restrains against ``ref.gro``, the structure as
built, exactly as the CHARMM-GUI README prescribes.
"""

import os
import shutil
import subprocess

MDP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mdp")

MINIMIZATION = "step6.0_minimization.mdp"
EQUILIBRATION = [f"step6.{i}_equilibration.mdp" for i in range(1, 7)]
PRODUCTION = "step7_production.mdp"


def stage_mdps(outdir):
    """Copy the packaged .mdp files into a run directory."""
    for name in [MINIMIZATION] + EQUILIBRATION + [PRODUCTION]:
        shutil.copy(os.path.join(MDP_DIR, name), os.path.join(outdir, name))


def find_gmx():
    """Locate the GROMACS executable, or explain what is missing."""
    exe = shutil.which("gmx") or shutil.which("gmx_mpi")
    if exe is None:
        raise RuntimeError(
            "no 'gmx' executable on PATH -- source your GMXRC "
            "(e.g. `source /usr/local/gromacs/bin/GMXRC`) and try again")
    return exe


def _stage(gmx, outdir, mdp, conf, out, cpt=None, ntomp=None, nsteps=None,
           maxwarn=5):
    grompp = [gmx, "grompp", "-f", mdp, "-c", conf, "-r", "ref.gro",
              "-p", "topol.top", "-n", "index.ndx", "-o", f"{out}.tpr",
              "-maxwarn", str(maxwarn)]
    if cpt:
        grompp += ["-t", cpt]
    subprocess.run(grompp, cwd=outdir, check=True)

    mdrun = [gmx, "mdrun", "-deffnm", out]
    if ntomp:
        mdrun += ["-ntmpi", "1", "-ntomp", str(ntomp)]
    if nsteps is not None:
        mdrun += ["-nsteps", str(nsteps)]
    subprocess.run(mdrun, cwd=outdir, check=True)
    return f"{out}.gro"


def run(outdir, production=False, ntomp=None, nsteps=None,
        structure="system.gro"):
    """Minimize, then equilibrate, a run directory staged by ``ld_lite build``.

    ``nsteps`` overrides the step count of every stage, which turns the whole
    protocol into a few-minute smoke test of the topology, index groups and
    restraints before committing to the real run.

    Returns the path to the final coordinate file.
    """
    for required in ("topol.top", "index.ndx", "ref.gro", structure):
        if not os.path.exists(os.path.join(outdir, required)):
            raise FileNotFoundError(
                f"{outdir}/{required} is missing -- run `ld_lite build` first")
    if not os.path.isdir(os.path.join(outdir, "toppar")):
        raise FileNotFoundError(
            f"{outdir}/toppar is missing -- copy your CHARMM36 toppar directory "
            f"there, or rebuild with `ld_lite build --toppar ...`")
    gmx = find_gmx()

    if nsteps is not None:
        print(f"[ld_lite] smoke test: every stage capped at {nsteps} steps")

    print("[ld_lite] step6.0  energy minimization")
    prev = _stage(gmx, outdir, MINIMIZATION, structure, "em",
                  ntomp=ntomp, nsteps=nsteps)

    # step6.1 generates its own velocities, so it starts without a checkpoint
    cpt = None
    for i, mdp in enumerate(EQUILIBRATION, start=1):
        print(f"[ld_lite] step6.{i}  equilibration")
        prev = _stage(gmx, outdir, mdp, prev, f"eq{i}", cpt=cpt,
                      ntomp=ntomp, nsteps=nsteps)
        cpt = f"eq{i}.cpt"

    if production:
        print("[ld_lite] step7   NPT production")
        prev = _stage(gmx, outdir, PRODUCTION, prev, "npt", cpt=cpt,
                      ntomp=ntomp, nsteps=nsteps)

    return os.path.join(outdir, prev)
