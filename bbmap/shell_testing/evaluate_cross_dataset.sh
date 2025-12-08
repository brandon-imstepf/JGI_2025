#!/bin/bash
#
# Evaluate July-trained networks on the October dataset (and validation subset).
# This script looks for .bbnet models produced by the July sweep (ncbi_sorted_*),
# then runs train.sh in evaluate mode against the October subset specified.
#
# Usage:
#   ./evaluate_cross_dataset.sh [subset] [max_models]
# Subset can be full|floats_only|onehot_only|floats_plus_onehot_no_middle (default: full).
# max_models (optional) limits how many models to evaluate (useful for a quick spot-check).
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_ROOT="$(cd "${BASE_DIR}/.." && pwd)"
TRAIN_SH="${BASE_DIR}/train.sh"

SUBSET="${1:-full}"
MAX_MODELS="${2:-0}"
ALLOWED_SUBSETS=(full floats_only onehot_only floats_plus_onehot_no_middle)
if [[ ! " ${ALLOWED_SUBSETS[*]} " =~ " ${SUBSET} " ]]; then
  echo "[ERROR] Invalid subset '${SUBSET}'. Choose one of: ${ALLOWED_SUBSETS[*]}" >&2
  exit 1
fi

JULY_MODELS_DIR="${PROJECT_ROOT}/comparison_models/hparam_sweep_nets/ncbi_sorted_${SUBSET}"
OCT_TRAIN="${PROJECT_ROOT}/training_datasets/NCBI_Datasets_Training_ncbi_truth_slurm_19888056_subsets/${SUBSET}/NCBI_Datasets_Training.ncbi_truth.all.subset_${SUBSET}.tsv.gz"
OCT_VALIDATE="${PROJECT_ROOT}/training_datasets/Validation_NCBI_10_ncbi_truth_20251120_122824_subsets/${SUBSET}/Validation_NCBI_10.ncbi_truth.all.subset_${SUBSET}.tsv.gz"

if [[ ! -d "${JULY_MODELS_DIR}" ]]; then
  echo "[ERROR] Model directory not found: ${JULY_MODELS_DIR}" >&2
  exit 1
fi
if [[ ! -f "${OCT_TRAIN}" ]]; then
  echo "[ERROR] October training subset missing: ${OCT_TRAIN}" >&2
  exit 1
fi
if [[ ! -f "${OCT_VALIDATE}" ]]; then
  echo "[ERROR] October validation subset missing: ${OCT_VALIDATE}" >&2
  exit 1
fi

RESULTS_CSV="${PROJECT_ROOT}/logs/hparam_sweep/cross_dataset_eval_${SUBSET}.csv"
echo "model,subset,ERR,WER,FPR,FNR,CTF" > "${RESULTS_CSV}"

count=0
for MODEL in "${JULY_MODELS_DIR}"/*.bbnet; do
  if [[ "${MAX_MODELS}" -gt 0 && "${count}" -ge "${MAX_MODELS}" ]]; then
    break
  fi
  MODEL_NAME="$(basename "${MODEL}")"
  LOG_PATH="${MODEL%.bbnet}_eval.log"
  echo "Evaluating ${MODEL_NAME} on October ${SUBSET} subset..."
  ${TRAIN_SH} evaluate=t netin="${MODEL}" in="${OCT_TRAIN}" validate="${OCT_VALIDATE}" >"${LOG_PATH}" 2>&1 || {
    echo "  [FAIL] See ${LOG_PATH}" >&2
    continue
  }
  BEST_LINE="$(grep '^Best:' "${LOG_PATH}" || true)"
  if [[ -z "${BEST_LINE}" ]]; then
    echo "  [WARN] No Best line found, see ${LOG_PATH}" >&2
    continue
  fi
  ERR=$(echo "${BEST_LINE}" | awk -F'ERR=' '{print $2}' | awk '{print $1}')
  WER=$(echo "${BEST_LINE}" | awk -F'WER=' '{print $2}' | awk '{print $1}')
  FPR=$(echo "${BEST_LINE}" | awk -F'FPR=' '{print $2}' | awk '{print $1}')
  FNR=$(echo "${BEST_LINE}" | awk -F'FNR=' '{print $2}' | awk '{print $1}')
  CTF=$(echo "${BEST_LINE}" | awk -F'CTF=' '{print $2}')
  echo "${MODEL_NAME},${SUBSET},${ERR},${WER},${FPR},${FNR},${CTF}" >> "${RESULTS_CSV}"
  count=$((count+1))
done

echo "Results saved to ${RESULTS_CSV}"
