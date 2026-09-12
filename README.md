# ld_lite

Build a **lipid droplet** starting structure from a phospholipid bilayer, then
run GROMACS minimization and equilibration on it. Two commands, pure Python for
the build step.

An LD is a neutral-lipid core wrapped by two phospholipid monolayers. `ld_lite`
makes one by taking a solvated, equilibrated CHARMM-GUI bilayer, prying the two
leaflets apart along z, and filling the vacated midplane with triolein (TRIO).
The bilayer's own water and ions ride along untouched.

```
    bilayer                          lipid droplet
  ~~~~~~~~~~~~  heads              ~~~~~~~~~~~~   monolayer
  ||||||||||||  tails      -->     ||||||||||||
  ||||||||||||  tails              ############   <- inserted TRIO core
  ~~~~~~~~~~~~  heads              ||||||||||||
                                   ~~~~~~~~~~~~   monolayer
```

## Install

```bash
pip install -e .          # or: pip install -e ".[test]"
```

Needs Python ≥3.8, `numpy`, `scipy`, `MDAnalysis`. GROMACS (`gmx` on `PATH`) is
needed only for the equilibration step, not for the build.

## Use

```bash
# 1. build the droplet and stage a run directory
python -m ld_lite build --bilayer bil.gro --n-trio 150 --out runs/ld1 \
    --toppar /path/to/charmm36/toppar

# 2. minimize + equilibrate it
python -m ld_lite equil runs/ld1
```

You bring two things: a **solvated phospholipid bilayer** `.gro` (CHARMM-GUI
output works directly — any composition, as long as the head groups have a `P`
atom) and your **CHARMM36 `toppar` directory**, which must contain an `.itp` for
every residue in the system including `TRIO.itp`.

Sanity-check the whole setup in a couple of minutes before committing to a real
run — this does every stage, just with the step count capped:

```bash
python -m ld_lite equil runs/ld1 --nsteps 50 --ntomp 8
```

`build` writes a self-contained run directory:

| file | |
|---|---|
| `system.gro` | the built droplet |
| `ref.gro` | copy of `system.gro`; every stage restrains against it |
| `topol.top` | `#include`s plus `[ molecules ]` counted in file order |
| `index.ndx` | `SOLU` (core) / `MEMB` (phospholipids + CHYO) / `SOLV` (water + ions) |
| `step6.*.mdp`, `step7_production.mdp` | the CHARMM-GUI protocol |
| `run_equil.sh`, `run_equil.slurm` | the same protocol in plain bash, for clusters |
| `toppar/` | if you passed `--toppar` |

`equil` then runs `step6.0` energy minimization followed by `step6.1`–`step6.6`,
which walk the position restraints from stiff down to free while the timestep
goes 1 → 2 fs and the barostat switches on. Add `--production` to continue into
unrestrained NPT (`step7`). Output is `em.gro`, `eq1.gro` … `eq6.gro`, `npt.gro`.

On a GPU node this is hours of wall time; `--ntomp N` sets the thread count.

## Running on a cluster (SLURM)

The Python `equil` command is a convenience for a local machine. For real runs
`build` also stages **`run_equil.sh`**, which does the identical protocol in
plain bash and needs *nothing but GROMACS* — no Python, no `ld_lite` install on
the cluster. Build locally, copy the directory over, submit:

```bash
rsync -a runs/ld1/ cluster:~/ld1/
ssh cluster
cd ~/ld1
sbatch run_equil.slurm          # or: bash run_equil.sh --production
```

`run_equil.slurm` is a **template** — the `--account`, `--partition` and
`--gres` lines are marked `FIXME` and must match your allocation. It sets the
knobs below and then calls `run_equil.sh`.

| variable | |
|---|---|
| `GMX` | GROMACS executable (default `gmx`) |
| `NT` | OpenMP threads (default `$SLURM_CPUS_PER_TASK`, else 8) |
| `MAXH` | stop cleanly after N hours and checkpoint — set just under `--time` |
| `MDRUN_OPTS` | extra mdrun flags, e.g. `-nb gpu -pme gpu` |
| `NSTEPS` | cap every stage, for a quick smoke test |

### It survives the wall clock

This is the part that matters on a queue. Every stage is **skipped if its
`.gro` already exists** and **resumed from its checkpoint if it was cut off**,
so a job that runs out of time is just resubmitted — nothing is recomputed:

```
$ sbatch run_equil.slurm          # hits --time during production
[ld_lite] npt: stopped before finishing (wall clock?).
[ld_lite] Resubmit this script to continue from npt.cpt.

$ sbatch run_equil.slurm          # picks up exactly where it stopped
[ld_lite] em:  already complete, skipping
[ld_lite] eq1: already complete, skipping
...
[ld_lite] npt: resuming from checkpoint
```

`MAXH` is what makes this work: mdrun stops itself and writes a checkpoint
*before* SLURM kills the job. Without it a hard kill can lose the whole stage.
The script exits 0 on a wall-clock stop (so the job counts as successful) but
propagates a real GROMACS failure as a non-zero exit.

## Cholesteryl ester in the leaflets

A bilayer that already contains **CHYO** carries straight through. Build the
CHYO-containing bilayer in CHARMM-GUI (which sizes the box for the real
per-lipid areas), then insert the core as usual:

```bash
python -m ld_lite build --bilayer chyo_bilayer.gro --n-trio 100 --out runs/ld_ce \
    --toppar /path/to/toppar
```

`ld_lite` only ever inserts TRIO, so any CHYO in the system is by definition a
**surface** lipid: it goes into the `MEMB` index group with the phospholipids,
not the `SOLU` core group, and the topology includes it as `toppar/chyo.itp`
(lowercase, which is how CHARMM-GUI ships it). Nothing to configure.

Do **not** try to reach a CHYO leaflet by substituting lipids in place — a bulky
substitution over-packs the leaflet and will not minimize. That is what
CHARMM-GUI Membrane Builder is for.

## Lipids it hasn't heard of

Composition comes entirely from your bilayer — `ld_lite` never builds or edits a
membrane, so any CHARMM-GUI mixture carries through. What it does need to know
is which thermostat group each residue belongs to. These are recognised out of
the box:

| group | residues |
|---|---|
| `MEMB` | POPC, DOPE, SAPI, POPE, POPG, POPS, PSM, CHL1, CHYO |
| `SOLV` | TIP3, SOL, SOD, CLA, POT, K, NA |
| `SOLU` | TRIO |

Anything else is declared on the command line — no source edit:

```bash
python -m ld_lite build --bilayer bil.gro --n-trio 150 --out runs/ld1 \
    --membrane DPPC DMPC POPI --solvent URE
```

Each declared residue joins that index group and is included as
`toppar/<RESNAME>.itp`. The same is available from Python as the `extra_core` /
`extra_membrane` / `extra_solvent` arguments to `write_index`, and `known=` to
`write_top`. An undeclared residue still fails loudly when the files are
written, naming the residue and the flag to use, rather than becoming a cryptic
`grompp` rejection later.

## Python API

```python
from ld_lite import build_ld, write_top, write_index

ld = build_ld("bil.gro", n_trio=150)
ld.atoms.write("system.gro")
write_top(ld, "topol.top")
write_index(ld, "index.ndx")
```

## How the core is placed

TRIO molecules are dropped onto a jittered 3-D grid inside the gap. Each one is
tried in 100 random orientations and the pose sitting furthest from everything
already placed wins — rigid-body relaxation cannot un-thread two interlocked
molecules, so the trick is to never interlock them. Afterwards the core is
relaxed as rigid bodies against the bilayer, which is held fixed: it is already
equilibrated, and pushing the lipids around only shoves them into water that
stays put.

`--packing-fraction` (default 0.25) is the main knob. It sets how tightly the
core is filled, and therefore how thick the inserted slab is. Lower is looser
and cleaner but leaves a taller box for NPT to compress back to real density.
Measured on a 90 POPC / 90 DOPE bilayer with 60 TRIO, counting inter-molecular
contacts closer than 1.2 Å in the final structure:

| `--packing-fraction` | contacts < 1.2 Å | closest contact |
|---|---|---|
| 0.35 | 119 | 0.39 Å |
| 0.25 (default) | 8 | 0.57 Å |
| 0.15 | 0 | 1.20 Å |

The build prints the closest contact it achieved. Anything above ~0.5 Å
minimizes without trouble under the restrained steepest descent in `step6.0`.

**Check that printed number rather than trusting the default**, because the
right fraction depends on how much floor area the box gives each molecule. The
table above comes from a roomy 118 × 118 Å box; packing the same 100 TRIO into
a 90 × 90 Å one, the default 0.25 gives only 0.44 Å and 0.15 is needed to reach
1.55 Å. A tighter box wants a lower fraction.

This is not a cosmetic number. A 60-TRIO build that came out at **0.46 Å**
minimized fine and cleared `step6.1`–`step6.3`, then died with LINCS constraint
errors at **`step6.4`** — the first stage at a 2 fs timestep — and did so
whether the stages were 25 steps or 200. Rebuilding the same system at
`--packing-fraction 0.15` gave **1.60 Å** and ran the whole protocol through
`step7` without a single warning. If the printed contact is much under ~1 Å,
lower the fraction and rebuild; the failure will otherwise surface several
stages later, long after minimization has convinced you it was fine.

## Limitations

- The output is a **starting structure**. Always minimize before running MD.
- Pure **TRIO** cores only. CHYO is supported in the *leaflets* (above), but a
  mixed triolein / cholesteryl-ester *core*, or a from-scratch grid build, needs
  the full [`ld_maker`](https://github.com/jaybraunjr/ld_maker) package.
- The bilayer must be solvated and ionized already — nothing here adds water.
- Head groups are found by atom name `P`, which covers PC/PE/PI/PS/PG.
- Molecules **split across the periodic boundary are reassembled** on load, and
  the count is printed. Trajectory frames written without `-pbc mol` routinely
  contain these; left alone they give a lipid a meaningless centre of mass (so
  it lands in the wrong leaflet) and make `grompp` abort with *"largest distance
  between excluded atoms … larger than the cut-off"* — an error `-maxwarn` will
  not absorb. Reassembly assumes a rectangular box and molecules shorter than
  half of it, which holds for CHARMM-GUI membranes.
- A residue with no `.itp` guess, or that lands in no index group, is reported
  when the files are written rather than as a cryptic `grompp` failure later.
  Extend `ITP_ORDER` / `MEMBRANE` in `ld_lite/gmx.py` for new lipids.
- `grompp` warns that `POSRES_FC_BB` / `POSRES_FC_SC` are unused. That is
  expected: those are protein backbone/sidechain restraints and these systems
  have no protein. `-maxwarn 5` absorbs it.

## Validation

Checked end to end with GROMACS 2022 on a CHARMM36 bilayer of 90 POPC / 90 DOPE
/ 31 CHYO plus water and ions (63,227 atoms), inserting a 100-TRIO core to give
a 79,927-atom droplet in a 118 × 118 × 127 Å box:

- closest inter-molecular contact in the built structure **1.46 Å**, with *no*
  contact under 1.0 Å; the core sits 100% between the leaflets and every CHYO
  stays in the leaflet slab
- all three thermostat groups populated — `SOLU` 39,699 DOF (TRIO core),
  `MEMB` 65,815 (phospholipids + CHYO), `SOLV` 71,421 (water + ions) — summing
  to every atom in the system
- `step6.0` steepest descents **converged to Fmax < 1000 in 221 of 5000 steps**,
  potential energy −5.83 × 10⁵ kJ/mol
- `step6.1`–`step6.6` ran with **zero LINCS warnings and no settle failures**.
  The box holds constant through the NVT stages and starts compressing once the
  semiisotropic barostat engages at `step6.3`, 1774.8 → 1749.8 nm³:

| | built | em | eq1 | eq2 | eq3 | eq4 | eq5 | eq6 |
|---|---|---|---|---|---|---|---|---|
| volume (nm³) | 1774.8 | 1774.8 | 1774.8 | 1774.8 | 1768.6 | 1758.9 | 1752.6 | 1749.8 |
| CHYO in leaflets | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |

Those equilibration stages were shortened to 1500 steps each to fit a CPU-only
test box — enough to show the protocol is stable, not enough to converge the
droplet. Run the full protocol for production.

## Tests

```bash
pytest -q
```

Synthetic fixtures only — fast, no large files, and no GROMACS required.
