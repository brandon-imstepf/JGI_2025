#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 12:00:00
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -J base_onehot
#SBATCH -o /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/hparam_baseline/baseline_onehot_only_%j.out
#SBATCH -e /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/logs/hparam_baseline/baseline_onehot_only_%j.err

set -euo pipefail

PROJECT_ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf"
TRAIN_SCRIPT="${PROJECT_ROOT}/bbmap/train.sh"
SUBSET="onehot_only"
INPUTS=348
MINDIMS="${INPUTS},128,96,64,32,1"
MAXDIMS="${INPUTS},256,192,128,64,1"
TRAIN_PATH="${PROJECT_ROOT}/training_datasets/NCBI_Datasets_10-23_SORTED_ncbi_truth_slurm_19888055_subsets/${SUBSET}/NCBI_Datasets_10-23_SORTED.ncbi_truth.all.subset_${SUBSET}.tsv.gz"
VALIDATE_PATH="${PROJECT_ROOT}/training_datasets/Validation_NCBI_10_ncbi_truth_20251120_122824_subsets/${SUBSET}/Validation_NCBI_10.ncbi_truth.all.subset_${SUBSET}.tsv.gz"
LOG_DIR="${PROJECT_ROOT}/logs/hparam_baseline"
NET_DIR="${PROJECT_ROOT}/comparison_models/hparam_baseline_nets"
LOG_PATH="${LOG_DIR}/baseline_${SUBSET}.log"
OUT_PATH="${NET_DIR}/baseline_${SUBSET}.bbnet"

mkdir -p "${LOG_DIR}" "${NET_DIR}"

${TRAIN_SCRIPT} -Xmx64g seed=12345 \
    in="${TRAIN_PATH}" \
    validate="${VALIDATE_PATH}" \
    mindims="${MINDIMS}" \
    maxdims="${MAXDIMS}" \
    nets=1 cycles=1 batches=200000 \
    out="${OUT_PATH}" \
    >"${LOG_PATH}" 2>&1
