#!/bin/bash
# Compile CATHY once (fills $PYCATHY_CACHE) before an ensemble: see submit_ensemble.sh
#SBATCH -J pycathy-build
#SBATCH -p generic              # [VERIFY] Drago asks to compile in the queue enabled for it: ask which one
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem 4G
#SBATCH -t 0-00:30
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err

set -euo pipefail
module purge
source "${PYCATHY_BASE:-$HOME/pycathy}/miniforge3/etc/profile.d/conda.sh"
conda activate pycathy
export PYCATHY_CACHE="${PYCATHY_CACHE:-$HOME/.cache/pycathy}"

WORK="${SLURM_TMPDIR:-/tmp}/pycathy_build_${SLURM_JOB_ID}"
python "${SLURM_SUBMIT_DIR}/run_case.py" "$WORK" --build-only   # [ADJUST] same model set-up => same CATHY.H
