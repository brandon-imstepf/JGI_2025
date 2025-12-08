#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
TRAIN_SCRIPT="${PROJECT_ROOT}/bbmap/train.sh"
LOG_ROOT="${PROJECT_ROOT}/logs/hparam_baseline"
NET_ROOT="${PROJECT_ROOT}/comparison_models/hparam_baseline_nets"

mkdir -p "${LOG_ROOT}" "${NET_ROOT}"

declare -A TRAIN_DIRS=(
  ["ncbi_sorted_ncbi"]="NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets"
  ["ncbi_sorted_ncbi_intersect"]="NCBI_Datasets_10-23_SORTED_intersect_truth_slurm_19888052_subsets"
  ["ncbi_training_ncbi"]="NCBI_Datasets_Training_ncbi_truth_slurm_19888056_subsets"
  ["ncbi_training_ncbi_intersect"]="NCBI_Datasets_Training_intersect_truth_slurm_19888054_subsets"
)

declare -A TRAIN_PREFIX=(
  ["ncbi_sorted_ncbi"]="NCBI_Datasets_10-23_SORTED.ncbi_truth"
  ["ncbi_sorted_ncbi_intersect"]="NCBI_Datasets_10-23_SORTED.intersect_truth"
  ["ncbi_training_ncbi"]="NCBI_Datasets_Training.ncbi_truth"
  ["ncbi_training_ncbi_intersect"]="NCBI_Datasets_Training.intersect_truth"
)

declare -A VALIDATE_DIRS=(
  ["ncbi"]="Validation_NCBI_10_ncbi_truth_20251120_122824_subsets"
  ["ncbi_intersect"]="Validation_NCBI_10_intersect_truth_20251120_123017_subsets"
)

declare -A VALIDATE_PREFIX=(
  ["ncbi"]="Validation_NCBI_10.ncbi_truth"
  ["ncbi_intersect"]="Validation_NCBI_10.intersect_truth"
)

declare -A INPUT_DIMS=(
  ["full"]=356
  ["floats_only"]=8
  ["onehot_only"]=348
  ["floats_plus_onehot_no_middle"]=240
)

loss_flags() {
  local name="$1"
  case "${name}" in
    mse)
      echo "loss=mse losslegacyweighting=t"
      ;;
    wbce)
      echo "loss=wbce losslegacyweighting=f lossposweight=1.0 lossnegweight=1.0"
      ;;
    tversky)
      echo "loss=tversky losslegacyweighting=f lossalpha=0.3 lossbeta=0.7 losssmooth=1e-6"
      ;;
    *)
      echo "Unknown loss type ${name}" >&2
      exit 1
      ;;
  esac
}

cohorts=(ncbi_sorted ncbi_training)
truths=(ncbi ncbi_intersect)
subsets=(full floats_only onehot_only floats_plus_onehot_no_middle)
losses=(mse wbce tversky)

for cohort in "${cohorts[@]}"; do
  for truth in "${truths[@]}"; do
    dataset_key="${cohort}_${truth}"
    train_dir="${TRAIN_DIRS[${dataset_key}]}"
    train_prefix="${TRAIN_PREFIX[${dataset_key}]}"
    validate_dir="${VALIDATE_DIRS[${truth}]}"
    validate_prefix="${VALIDATE_PREFIX[${truth}]}"

    if [[ -z "${train_dir:-}" || -z "${validate_dir:-}" ]]; then
      echo "Missing directory mapping for dataset=${cohort}, truth=${truth}" >&2
      exit 1
    fi

    for subset in "${subsets[@]}"; do
      input_dim="${INPUT_DIMS[${subset}]}"
      if [[ -z "${input_dim:-}" ]]; then
        echo "Unknown subset ${subset}" >&2
        exit 1
      fi

      mindims="${input_dim},128,96,64,32,1"
      maxdims="${input_dim},256,192,128,64,1"

      train_path="${PROJECT_ROOT}/training_datasets/${train_dir}/${subset}/${train_prefix}.all.subset_${subset}.tsv.gz"
      validate_path="${PROJECT_ROOT}/training_datasets/${validate_dir}/${subset}/${validate_prefix}.all.subset_${subset}.tsv.gz"

      if [[ ! -f "${train_path}" ]]; then
        echo "Missing training file: ${train_path}" >&2
        exit 1
      fi
      if [[ ! -f "${validate_path}" ]]; then
        echo "Missing validation file: ${validate_path}" >&2
        exit 1
      fi

      combo_dir="${LOG_ROOT}/${cohort}_${truth}_${subset}"
      net_dir="${NET_ROOT}/${cohort}_${truth}_${subset}"
      mkdir -p "${combo_dir}" "${net_dir}"

      for loss in "${losses[@]}"; do
        loss_args="$(loss_flags "${loss}")"
        job_label="baseline_${cohort}_${truth}_${subset}_${loss}"
        log_file="${combo_dir}/${job_label}.log"
        net_path="${net_dir}/${job_label}.bbnet"
        slurm_out="${combo_dir}/${job_label}.slurm.out"
        slurm_err="${combo_dir}/${job_label}.slurm.err"
        wrap_cmd="${TRAIN_SCRIPT} -Xmx64g seed=12345 in=${train_path} validate=${validate_path} mindims=${mindims} maxdims=${maxdims} nets=1 cycles=1 batches=200000 out=${net_path} ${loss_args} > ${log_file} 2>&1"

        echo "Submitting ${job_label}"
        sbatch -A grp-org-gt-ga -q jgi_normal -t 12:00:00 -c 16 --mem=64G \
          --job-name="base_${subset}" \
          --output="${slurm_out}" \
          --error="${slurm_err}" \
          --wrap="${wrap_cmd}"
        sleep 1
      done
    done
  done
done
