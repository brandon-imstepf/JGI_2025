#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 72:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=ncbi_truth_training
#SBATCH --output=ncbi_truth_training_%j.log
#SBATCH --error=ncbi_truth_training_%j.err

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
SHELL_TESTING_DIR="${PROJECT_ROOT}/bbmap/shell_testing"
DATASET_DIR="${PROJECT_ROOT}/NCBI_Datasets_Training"

echo "[$(date)] Starting NCBI-truth dataset build for ${DATASET_DIR}"

cd "${SHELL_TESTING_DIR}"

export RUN_TAG="${RUN_TAG:-slurm_${SLURM_JOB_ID:-manual}}"

./build_ncbi_truth_dataset.sh \
    in="${DATASET_DIR}" \
    label="NCBI_Datasets_Training" \
    run_tag="${RUN_TAG}"

echo "[$(date)] Completed NCBI-truth dataset build for NCBI_Datasets_Training"
