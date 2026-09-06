"""Command line: build a lipid droplet, then minimize and equilibrate it.

    python -m ld_lite build --bilayer bil.gro --n-trio 150 --out runs/ld1 \
        --toppar /path/to/toppar
    python -m ld_lite equil runs/ld1
"""

import argparse
import os
import shutil
import subprocess

from .build import build_ld
from .gmx import write_top, write_index
from . import equil


def _build(args):
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)

    print(f"[ld_lite] inserting {args.n_trio} TRIO into {args.bilayer}")
    u = build_ld(args.bilayer, args.n_trio,
                 core_source=args.core_source,
                 core_thickness=args.core_thickness,
                 packing_fraction=args.packing_fraction,
                 seed=args.seed,
                 relax=not args.no_relax,
                 d_min=args.d_min)

    gro = os.path.join(outdir, "system.gro")
    u.atoms.write(gro)
    declared = args.membrane + args.solvent + args.core
    write_top(u, os.path.join(outdir, "topol.top"), known=declared)
    counts = write_index(u, os.path.join(outdir, "index.ndx"),
                         extra_membrane=args.membrane,
                         extra_solvent=args.solvent,
                         extra_core=args.core)
    # the restraints in every stage reference the structure as built
    shutil.copy(gro, os.path.join(outdir, "ref.gro"))
    equil.stage_mdps(outdir)

    if args.toppar:
        dest = os.path.join(outdir, "toppar")
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(args.toppar, dest)

    print(f"[ld_lite] wrote {gro}  ({u.atoms.n_atoms} atoms)")
    print(f"[ld_lite] box = {u.dimensions[:3].round(1)} A")
    print("[ld_lite] index groups: " +
          ", ".join(f"{k}={v}" for k, v in counts.items()))
    if not args.toppar:
        print(f"[ld_lite] note: copy your CHARMM36 toppar/ into {outdir} before running gmx")
    print(f"[ld_lite] next: python -m ld_lite equil {outdir}")


def _equil(args):
    final = equil.run(args.rundir, production=args.production,
                      ntomp=args.ntomp, nsteps=args.nsteps)
    print(f"[ld_lite] done -> {final}")


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="ld_lite",
        description="Build a lipid droplet from a bilayer, then minimize and equilibrate it.")
    sub = p.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="insert a TRIO core into a bilayer and stage a run directory")
    b.add_argument("--bilayer", required=True, help="solvated phospholipid bilayer .gro")
    b.add_argument("--n-trio", type=int, required=True, help="TRIO molecules in the core")
    b.add_argument("--out", required=True, help="run directory to create")
    b.add_argument("--toppar", default=None, help="CHARMM36 toppar directory to copy in")
    b.add_argument("--core-source", default=None,
                   help="structure to take the TRIO conformation from (default: packaged)")
    b.add_argument("--core-thickness", type=float, default=None,
                   help="core slab thickness in A (default: from --packing-fraction)")
    b.add_argument("--packing-fraction", type=float, default=0.25,
                   help="initial core fill, 0-1; lower leaves fewer overlaps "
                        "but a taller box for NPT to compress (default 0.25)")
    b.add_argument("--d-min", type=float, default=1.6,
                   help="target minimum inter-molecular contact, A")
    b.add_argument("--no-relax", action="store_true", help="skip clash relaxation")
    b.add_argument("--seed", type=int, default=0, help="RNG seed")
    b.add_argument("--membrane", nargs="+", default=[], metavar="RESNAME",
                   help="extra residues to treat as membrane lipids, e.g. "
                        "--membrane DPPC DMPC. Each is included as "
                        "toppar/<RESNAME>.itp and joins the MEMB group")
    b.add_argument("--solvent", nargs="+", default=[], metavar="RESNAME",
                   help="extra residues to treat as solvent or ions (SOLV group)")
    b.add_argument("--core", nargs="+", default=[], metavar="RESNAME",
                   help="extra residues to treat as core neutral lipids (SOLU group)")
    b.set_defaults(func=_build)

    e = sub.add_parser("equil", help="run step6.0 minimization + step6.1-6.6 equilibration")
    e.add_argument("rundir", help="run directory created by `ld_lite build`")
    e.add_argument("--production", action="store_true", help="also run step7 NPT production")
    e.add_argument("--ntomp", type=int, default=None, help="OpenMP threads for mdrun")
    e.add_argument("--nsteps", type=int, default=None,
                   help="cap every stage at N steps -- a fast smoke test of the "
                        "topology and index groups before the real run")
    e.set_defaults(func=_equil)

    args = p.parse_args(argv)
    try:
        args.func(args)
    except (RuntimeError, FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"[ld_lite] error: {exc}")
    except subprocess.CalledProcessError as exc:
        raise SystemExit(
            f"[ld_lite] error: {' '.join(exc.cmd[:2])} failed with exit code "
            f"{exc.returncode} -- see the GROMACS output above")


if __name__ == "__main__":
    main()
