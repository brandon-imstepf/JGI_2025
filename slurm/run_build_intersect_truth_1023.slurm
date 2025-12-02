#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 72:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=intersect_truth_1023
#SBATCH --output=intersect_truth_1023_%j.log
#SBATCH --error=intersect_truth_1023_%j.err

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
SHELL_TESTING_DIR="${PROJECT_ROOT}/bbmap/shell_testing"
DATASET_DIR="${PROJECT_ROOT}/NCBI_Datasets_10-23_SORTED"

echo "[$(date)] Starting intersect-truth dataset build for ${DATASET_DIR}"

cd "${SHELL_TESTING_DIR}"

export RUN_TAG="${RUN_TAG:-slurm_${SLURM_JOB_ID:-manual}}"

./build_intersect_truth_dataset.sh \
    in="${DATASET_DIR}" \
    label="NCBI_Datasets_10-23_SORTED" \
    run_tag="${RUN_TAG}"

echo "[$(date)] Completed intersect-truth dataset build for NCBI_Datasets_10-23_SORTED"
