#!/bin/bash
# One pyCATHY run (e.g. your monthly LAI loop). Submit from the folder holding run_case.py:
#   sbatch pycathy_single.sbatch
#SBATCH -J pycathy
#SBATCH -p generic              # CPU partition (gpu / special also exist). [VERIFY] limits with the Drago team
#SBATCH -N 1
#SBATCH -n 1                    # CATHY is serial: one core
#SBATCH --mem 8G                # [ADJUST] a bare number means MB; check real use with `seff <jobid>`
#SBATCH -t 0-12:00              # D-HH:MM
#SBATCH -o slurm.%x.%j.out
#SBATCH -e slurm.%x.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@csic.es

set -euo pipefail
module purge                    # [VERIFY] no module needed: compiler and Python come from conda
source "${PYCATHY_BASE:-$HOME/pycathy}/miniforge3/etc/profile.d/conda.sh"   # [ADJUST] if conda comes from a module
conda activate pycathy

# Build cache shared by all your jobs (Lustre group/home dir): compiled once per CATHY.H.
export PYCATHY_CACHE="${PYCATHY_CACHE:-$HOME/.cache/pycathy}"

# CATHY writes many output files: run on the node-local job directory, copy back at the end.
WORK="${SLURM_TMPDIR:-/tmp}/pycathy_${SLURM_JOB_ID}"
RESULTS="${SLURM_SUBMIT_DIR}/results_${SLURM_JOB_ID}"
mkdir -p "$WORK" "$RESULTS"
trap 'rsync -a "$WORK"/ "$RESULTS"/' EXIT      # also on failure / SIGTERM at the time limit

python "${SLURM_SUBMIT_DIR}/run_case.py" "$WORK"
