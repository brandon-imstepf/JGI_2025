#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 72:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=nn_sensitivity
#SBATCH --output=nn_sensitivity_%j.log
#SBATCH --error=nn_sensitivity_%j.err

# Optional environment overrides before sbatch:
#   export NN_STRENGTH_VALUES="0 0.5 1"
#   export NN_CUTOFF_VALUES="0.3 0.4 0.5 0.6"
#   export NN_MULT_MODE=ratio_pow2
#   export NN_PREFILTER=0.0
#   export GENOME_LIMIT=3

set -euo pipefail

echo "[$(date)] Starting NN sensitivity sweep job (SLURM_JOB_ID=${SLURM_JOB_ID:-manual})"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHELL_TESTING_DIR="$(cd "${SCRIPT_DIR}/../bbmap/shell_testing" && pwd)"
cd "${SHELL_TESTING_DIR}"

export RUN_TAG="${RUN_TAG:-slurm_${SLURM_JOB_ID:-manual}}"
export NN_STRENGTH_VALUES="${NN_STRENGTH_VALUES:-"0 0.5 1.0 2.0 5.0 10.0"}"
export NN_CUTOFF_VALUES="${NN_CUTOFF_VALUES:-"0.3 0.4 0.5 0.6 0.7 0.8 0.9"}"
export NN_MULT_MODE="${NN_MULT_MODE:-ratio}"
export NN_PREFILTER="${NN_PREFILTER:-0.0}"
export GENOME_LIMIT="${GENOME_LIMIT:-}"

./run_nn_sensitivity.sh

echo "[$(date)] Finished NN sensitivity sweep job"
