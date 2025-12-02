#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=feature_subsets
#SBATCH --output=feature_subsets_%j.log
#SBATCH --error=feature_subsets_%j.err

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
SHELL_TESTING_DIR="${PROJECT_ROOT}/bbmap/shell_testing"
SCRIPT="${SHELL_TESTING_DIR}/build_feature_subsets.py"

declare -a DATASETS=(
    "NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055"
    "NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052"
    "NCBI_Datasets_Training_ncbi_truth_slurm_19888056"
    "NCBI_Datasets_Training_intersect_truth_slurm_19888054"
)

echo "[$(date)] Starting feature subset generation (SLURM_JOB_ID=${SLURM_JOB_ID:-manual})"

cd "${SHELL_TESTING_DIR}"

for dataset in "${DATASETS[@]}"; do
    IN_DIR="${PROJECT_ROOT}/training_datasets/${dataset}"
    OUT_DIR="${PROJECT_ROOT}/training_datasets/${dataset}_subsets"
    echo "[$(date)] Building subsets for ${dataset}"
    python3 "${SCRIPT}" --input-dir "${IN_DIR}" --output-dir "${OUT_DIR}"
done

echo "[$(date)] Finished feature subset generation"
