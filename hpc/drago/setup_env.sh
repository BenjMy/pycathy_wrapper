#!/usr/bin/env bash
# One-time install of pyCATHY on Drago (CSIC). Run it inside an interactive job, NOT on the login node:
#   srun -p generic -N1 -n1 --mem=8G -t 0-01:00 --pty bash -i
#   bash setup_env.sh
# Untested on Drago: check the [VERIFY] points of README.md first.
set -euo pipefail

BASE="${PYCATHY_BASE:-$HOME/pycathy}"          # Lustre home: /lustre/home/<centro>/<usuario>/pycathy
BRANCH="${PYCATHY_BRANCH:-refactor/no-gfortran-install}"
mkdir -p "$BASE"

if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"                # conda provided by a module or already installed
else
    # [VERIFY] `module spider conda` / `module spider miniforge` first: prefer the site module.
    ARCH="$(uname -m)"                             # x86_64 or aarch64
    curl -fL -o "$BASE/miniforge.sh" \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh"
    bash "$BASE/miniforge.sh" -b -p "$BASE/miniforge3"
    source "$BASE/miniforge3/etc/profile.d/conda.sh"
fi

conda create -y -n pycathy python=3.10
conda activate pycathy
conda install -y -c conda-forge fortran-compiler   # gfortran from conda: no module, no BLAS/LAPACK needed
pip install "git+https://github.com/BenjMy/pycathy_wrapper@${BRANCH}"

python - <<'PY'
from pyCATHY import fortran_build as fb
print("Fortran compiler:", fb.find_fortran_compiler())
print("flags:", fb.default_compile_flags(compiler=fb.find_fortran_compiler()))
PY
echo "OK. Quota: lfs quota -h \$HOME"
