#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 72:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=build_training_v2
#SBATCH --output=build_training_v2_%j.log
#SBATCH --error=build_training_v2_%j.err

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
SHELL_TESTING_DIR="${PROJECT_ROOT}/bbmap/shell_testing"
OUT_ROOT="${PROJECT_ROOT}/training_datasets_v2"

mkdir -p "${OUT_ROOT}"

declare -a DATASET_NAMES=(
    "NCBI_Datasets_Training"
    "NCBI_Datasets_10-23_SORTED"
)

echo "[$(date)] Starting dataset rebuild into ${OUT_ROOT}"
cd "${SHELL_TESTING_DIR}"

# Allow callers to override RUN_TAG for all invocations; otherwise stamp once.
GLOBAL_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"

for dataset in "${DATASET_NAMES[@]}"; do
    INPUT_DIR="${PROJECT_ROOT}/${dataset}"
    if [[ ! -d "${INPUT_DIR}" ]]; then
        echo "[$(date)] [WARN] Input directory ${INPUT_DIR} not found, skipping"
        continue
    fi

    echo "[$(date)] Building ncbi_truth dataset for ${dataset}"
    ./build_ncbi_truth_dataset.sh \
        in="${INPUT_DIR}" \
        outdir="${OUT_ROOT}" \
        label="${dataset}" \
        run_tag="${GLOBAL_TAG}_${dataset}_ncbi_truth"

    echo "[$(date)] Building intersect_truth dataset for ${dataset}"
    ./build_intersect_truth_dataset.sh \
        in="${INPUT_DIR}" \
        outdir="${OUT_ROOT}" \
        label="${dataset}" \
        run_tag="${GLOBAL_TAG}_${dataset}_intersect_truth"
done

echo "[$(date)] Dataset rebuild complete."
