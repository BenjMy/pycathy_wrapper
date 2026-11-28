#!/bin/bash
# DA ensemble / parameter sweep: one array task = one member. See submit_ensemble.sh
#SBATCH -J pycathy-ens
#SBATCH -p generic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8G                # [ADJUST]
#SBATCH -t 0-12:00
#SBATCH -o slurm.%x.%A_%a.out
#SBATCH -e slurm.%x.%A_%a.err

set -euo pipefail
module purge
source "${PYCATHY_BASE:-$HOME/pycathy}/miniforge3/etc/profile.d/conda.sh"
conda activate pycathy
export PYCATHY_CACHE="${PYCATHY_CACHE:-$HOME/.cache/pycathy}"   # atomic writes: safe with concurrent tasks

MEMBER="${SLURM_ARRAY_TASK_ID}"
WORK="${SLURM_TMPDIR:-/tmp}/pycathy_${SLURM_ARRAY_JOB_ID}_${MEMBER}"
RESULTS="${SLURM_SUBMIT_DIR}/results_${SLURM_ARRAY_JOB_ID}/member_$(printf '%03d' "$MEMBER")"
mkdir -p "$WORK" "$RESULTS"
trap 'rsync -a "$WORK"/ "$RESULTS"/' EXIT

python "${SLURM_SUBMIT_DIR}/run_case.py" "$WORK" --member "$MEMBER"
