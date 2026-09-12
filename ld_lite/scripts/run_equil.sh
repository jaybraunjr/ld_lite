#!/bin/bash
# =============================================================================
# ld_lite: CHARMM-GUI minimization + equilibration, in plain bash.
#
# Self-contained -- needs only GROMACS. No Python, no ld_lite install. Copy the
# whole run directory to the cluster and run this from inside it.
#
#   bash run_equil.sh                 # step6.0 -> step6.6
#   bash run_equil.sh --production    # ... and on into step7 NPT
#
# Environment knobs (all optional):
#   GMX=gmx_mpi        GROMACS executable                  (default: gmx)
#   NT=16              OpenMP threads   (default: $SLURM_CPUS_PER_TASK, else 8)
#   MAXH=3.9           stop cleanly after N hours of wall clock, writing a
#                      checkpoint. Set this just under your SLURM --time and
#                      resubmit to continue; see "Resuming" below.
#   MDRUN_OPTS="-nb gpu -pme gpu"      extra mdrun flags, e.g. GPU offload
#   NSTEPS=50          cap every stage at N steps -- a fast smoke test that the
#                      topology and index groups are right
#
# Resuming: every stage is skipped if its .gro already exists, and continued
# from its checkpoint if it was cut off part way. So if the job hits the wall
# clock, just submit it again -- it picks up exactly where it stopped.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

GMX=${GMX:-gmx}
NT=${NT:-${SLURM_CPUS_PER_TASK:-8}}
MAXH=${MAXH:-}
MDRUN_OPTS=${MDRUN_OPTS:-}
NSTEPS=${NSTEPS:-}
PRODUCTION=0
[ "${1:-}" = "--production" ] && PRODUCTION=1

command -v "$GMX" >/dev/null 2>&1 || {
    echo "ERROR: '$GMX' not found. Load your GROMACS module or source GMXRC," >&2
    echo "       or set GMX=/path/to/gmx." >&2
    exit 1; }

for f in topol.top index.ndx ref.gro system.gro; do
    [ -f "$f" ] || { echo "ERROR: $f missing -- run 'ld_lite build' first." >&2; exit 1; }
done
[ -d toppar ] || { echo "ERROR: toppar/ missing -- copy your CHARMM36 toppar here." >&2; exit 1; }

mdrun_flags=( -ntmpi 1 -ntomp "$NT" )
[ -n "$MAXH" ]   && mdrun_flags+=( -maxh "$MAXH" )
[ -n "$NSTEPS" ] && mdrun_flags+=( -nsteps "$NSTEPS" )
# shellcheck disable=SC2206
[ -n "$MDRUN_OPTS" ] && mdrun_flags+=( $MDRUN_OPTS )

# run_stage <name> <mdp> <input.gro> [input.cpt]
run_stage() {
    local name=$1 mdp=$2 conf=$3 cpt=${4:-}

    if [ -f "${name}.gro" ]; then
        echo "[ld_lite] ${name}: already complete, skipping"
        return 0
    fi

    if [ -f "${name}.cpt" ]; then
        echo "[ld_lite] ${name}: resuming from checkpoint ($(date))"
        "$GMX" mdrun -deffnm "$name" -cpi "${name}.cpt" "${mdrun_flags[@]}"
    else
        echo "[ld_lite] ${name}: $mdp ($(date))"
        local gro=( -f "$mdp" -c "$conf" -r ref.gro -p topol.top -n index.ndx
                    -o "${name}.tpr" -maxwarn 5 )
        # -maxwarn 5 absorbs the expected POSRES_FC_BB/SC "macro unused" note:
        # those are protein restraints and these systems have no protein.
        [ -n "$cpt" ] && [ -f "$cpt" ] && gro+=( -t "$cpt" )
        "$GMX" grompp "${gro[@]}"
        "$GMX" mdrun -deffnm "$name" "${mdrun_flags[@]}"
    fi

    if [ ! -f "${name}.gro" ]; then
        echo "[ld_lite] ${name}: stopped before finishing (wall clock?)."
        echo "[ld_lite] Resubmit this script to continue from ${name}.cpt."
        exit 0
    fi
}

echo "=== ld_lite equilibration | $GMX | ${NT} threads | $(date) ==="
[ -n "$NSTEPS" ] && echo "=== SMOKE TEST: every stage capped at ${NSTEPS} steps ==="

# step6.0 minimization, from the structure as built
run_stage em "step6.0_minimization.mdp" system.gro

# step6.1 generates its own velocities, so it starts without a checkpoint;
# every later stage continues from the previous one's checkpoint.
prev=em
cpt=""
for i in 1 2 3 4 5 6; do
    run_stage "eq${i}" "step6.${i}_equilibration.mdp" "${prev}.gro" "$cpt"
    prev="eq${i}"
    cpt="eq${i}.cpt"
done

if [ "$PRODUCTION" = "1" ]; then
    run_stage npt "step7_production.mdp" "${prev}.gro" "$cpt"
    prev=npt
fi

echo "=== COMPLETE -> ${prev}.gro | $(date) ==="
