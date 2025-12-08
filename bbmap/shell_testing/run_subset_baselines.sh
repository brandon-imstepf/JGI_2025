#!/bin/bash
# Run default training configuration on each feature subset (full, floats_only, onehot_only, no_middle)
# using the July ncbi-truth dataset. Produces baseline .bbnet files and logs
# so we can pick the best-performing subset before sweeping parameters.

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
TRAIN_SCRIPT="${PROJECT_ROOT}/bbmap/train.sh"
TRAIN_BASE="${PROJECT_ROOT}/training_datasets/NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets"
VALIDATE_BASE="${PROJECT_ROOT}/training_datasets/Validation_NCBI_10_ncbi_truth_20251120_122824_subsets"
TRAIN_PREFIX="NCBI_Datasets_10-23_SORTED.ncbi_truth"
VALIDATE_PREFIX="Validation_NCBI_10.ncbi_truth"
LOG_DIR="${PROJECT_ROOT}/logs/hparam_baseline"
NET_DIR="${PROJECT_ROOT}/comparison_models/hparam_baseline_nets"
SUBSETS=(full floats_only onehot_only floats_plus_onehot_no_middle)

mkdir -p "${LOG_DIR}" "${NET_DIR}"

function input_dimension() {
  case "$1" in
    full) echo 356 ;;
    floats_only) echo 8 ;;
    onehot_only) echo 348 ;;
    floats_plus_onehot_no_middle) echo 240 ;;
    *) echo "Unknown subset $1" >&2; exit 1 ;;
  esac
}

for subset in "${SUBSETS[@]}"; do
  train_file="${TRAIN_BASE}/${subset}/${TRAIN_PREFIX}.all.subset_${subset}.tsv.gz"
  validate_file="${VALIDATE_BASE}/${subset}/${VALIDATE_PREFIX}.all.subset_${subset}.tsv.gz"
  if [[ ! -f "${train_file}" ]]; then
    echo "[WARN] missing train file for ${subset}: ${train_file}" >&2
    continue
  fi
  if [[ ! -f "${validate_file}" ]]; then
    echo "[WARN] missing validate file for ${subset}: ${validate_file}" >&2
    continue
  fi
  inputs=$(input_dimension "${subset}")
  mindims="${inputs},128,96,64,32,1"
  maxdims="${inputs},256,192,128,64,1"
  log_name="baseline_${subset}"
  log_path="${LOG_DIR}/${log_name}.log"
  out_path="${NET_DIR}/${log_name}.bbnet"
  echo "Running baseline for subset=${subset} (inputs=${inputs})"
  "${TRAIN_SCRIPT}" -Xmx64g seed=12345 \
    in="${train_file}" \
    validate="${validate_file}" \
    mindims="${mindims}" \
    maxdims="${maxdims}" \
    nets=1 cycles=1 batches=200000 \
    out="${out_path}" \
    >"${log_path}" 2>&1
  echo "  -> log: ${log_path}"
  echo "  -> net: ${out_path}"
done
