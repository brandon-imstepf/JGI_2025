#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 12:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -J base_subset
#SBATCH -o /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/hparam_baseline/base_%x_%j.out
#SBATCH -e /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/hparam_baseline/base_%x_%j.err

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: sbatch baseline_subset.slurm <dataset_name> <subset> <truth_mode>" >&2
  echo " dataset_name: ncbi_sorted | ncbi_sorted_intersect | ncbi_training | ncbi_training_intersect" >&2
  echo " subset: full | floats_only | onehot_only | floats_plus_onehot_no_middle" >&2
  echo " truth_mode: ncbi | ncbi_intersect" >&2
  exit 1
fi

DATASET_NAME="$1"
SUBSET="$2"
TRUTH="$3"

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
TRAIN_SCRIPT="${PROJECT_ROOT}/bbmap/train.sh"

declare -A TRAIN_DIRS=(
  ["ncbi_sorted"]="NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets"
  ["ncbi_sorted_intersect"]="NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052_subsets"
  ["ncbi_training"]="NCBI_Datasets_Training_ncbi_truth_slurm_19888056_subsets"
  ["ncbi_training_intersect"]="NCBI_Datasets_Training_intersect_truth_slurm_19888054_subsets"
)

declare -A TRAIN_PREFIX=(
  ["ncbi_sorted"]="NCBI_Datasets_10-23_SORTED.ncbi_truth"
  ["ncbi_sorted_intersect"]="NCBI_Datasets_10-23_SORTED.intersect_truth"
  ["ncbi_training"]="NCBI_Datasets_Training.ncbi_truth"
  ["ncbi_training_intersect"]="NCBI_Datasets_Training.intersect_truth"
)

declare -A INPUTS=(
  ["full"]=356
  ["floats_only"]=8
  ["onehot_only"]=348
  ["floats_plus_onehot_no_middle"]=240
)

declare -A VALIDATE_DIRS=(
  ["ncbi"]="Validation_NCBI_10_ncbi_truth_20251120_122824_subsets"
  ["ncbi_intersect"]="Validation_NCBI_10_intersect_truth_20251120_123017_subsets"
)

declare -A VALIDATE_PREFIX=(
  ["ncbi"]="Validation_NCBI_10.ncbi_truth"
  ["ncbi_intersect"]="Validation_NCBI_10.intersect_truth"
)

if [[ -z "${TRAIN_DIRS[${DATASET_NAME}]:-}" ]]; then
  echo "Unknown dataset ${DATASET_NAME}" >&2
  exit 1
fi
if [[ -z "${INPUTS[${SUBSET}]:-}" ]]; then
  echo "Unknown subset ${SUBSET}" >&2
  exit 1
fi
if [[ -z "${VALIDATE_DIRS[${TRUTH}]:-}" ]]; then
  echo "Unknown truth mode ${TRUTH}" >&2
  exit 1
fi

TRAIN_PATH="${PROJECT_ROOT}/training_datasets/${TRAIN_DIRS[${DATASET_NAME}]}/${SUBSET}/${TRAIN_PREFIX[${DATASET_NAME}]}.all.subset_${SUBSET}.tsv.gz"
VALIDATE_PATH="${PROJECT_ROOT}/training_datasets/${VALIDATE_DIRS[${TRUTH}]}/${SUBSET}/${VALIDATE_PREFIX[${TRUTH}]}.all.subset_${SUBSET}.tsv.gz"

if [[ ! -f "${TRAIN_PATH}" ]]; then
  echo "Missing training file: ${TRAIN_PATH}" >&2
  exit 1
fi
if [[ ! -f "${VALIDATE_PATH}" ]]; then
  echo "Missing validation file: ${VALIDATE_PATH}" >&2
  exit 1
fi

INPUT_DIM=${INPUTS[${SUBSET}]}
MINDIMS="${INPUT_DIM},128,96,64,32,1"
MAXDIMS="${INPUT_DIM},256,192,128,64,1"

LOG_DIR="${PROJECT_ROOT}/logs/hparam_baseline"
NET_DIR="${PROJECT_ROOT}/comparison_models/hparam_baseline_nets"
BASE_NAME="baseline_${DATASET_NAME}_${SUBSET}_${TRUTH}"
LOG_PATH="${LOG_DIR}/${BASE_NAME}.log"
OUT_PATH="${NET_DIR}/${BASE_NAME}.bbnet"

mkdir -p "${LOG_DIR}" "${NET_DIR}"

echo "Running baseline for dataset=${DATASET_NAME}, subset=${SUBSET}, truth=${TRUTH}"

"${TRAIN_SCRIPT}" -Xmx64g seed=12345 \
    in="${TRAIN_PATH}" \
    validate="${VALIDATE_PATH}" \
    mindims="${MINDIMS}" \
    maxdims="${MAXDIMS}" \
    nets=1 cycles=1 batches=200000 \
    out="${OUT_PATH}" \
    > "${LOG_PATH}" 2>&1
