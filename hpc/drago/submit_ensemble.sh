#!/usr/bin/env bash
# Usage: bash submit_ensemble.sh 50 10      # 50 members, at most 10 running at once
set -euo pipefail
N="${1:?number of members}"; MAXRUN="${2:-10}"
BUILD=$(sbatch --parsable pycathy_warmcache.sbatch)
sbatch --dependency=afterok:"$BUILD" --array=0-$((N-1))%"$MAXRUN" pycathy_array.sbatch
