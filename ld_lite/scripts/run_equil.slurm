#!/bin/bash
# =============================================================================
# ld_lite equilibration on SLURM.   sbatch run_equil.slurm
#
# This is a TEMPLATE -- the account, partition and GPU type below are
# site-specific. Check them against your own allocation before submitting.
#
# If the job hits the wall clock it stops cleanly and checkpoints (see MAXH
# below), so simply `sbatch run_equil.slurm` again to continue where it left
# off. Nothing is recomputed.
# =============================================================================
#SBATCH --job-name=ld_equil
#SBATCH --output=ld_equil-%j.out
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1                 # FIXME: e.g. gpu:2080ti:1
#SBATCH --account=CHANGE-ME          # FIXME: e.g. swanson-gpu-np
#SBATCH --partition=CHANGE-ME        # FIXME: e.g. swanson-gpu-np

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- GROMACS ---------------------------------------------------------------
# FIXME: use whatever your cluster provides. `module spider gromacs` to find it.
module load gcc/8.5.0 openmpi/4.1.1 || true
module load gromacs               || true
# If GROMACS is not a module, source your build instead:
# source /path/to/gromacs/bin/GMXRC

# --- knobs passed through to run_equil.sh ----------------------------------
export NT=${SLURM_CPUS_PER_TASK:-8}

# Stop ~10 min before SLURM kills the job, so mdrun writes a checkpoint and the
# next submission resumes instead of losing the stage. Keep this under --time.
export MAXH=3.8

# GPU offload. Drop this line for a CPU-only partition.
export MDRUN_OPTS="-nb gpu -pme gpu"

bash run_equil.sh "$@"
