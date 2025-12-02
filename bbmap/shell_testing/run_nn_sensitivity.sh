#!/bin/bash
set -euo pipefail
trap 'echo "[ERROR] Command failed on line ${LINENO}. Exiting." >&2' ERR

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BBMAP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CURRENT_DIR="${BBMAP_DIR}/current"
BASE_DIR="$(cd "${BBMAP_DIR}/.." && pwd)"
INPUT_DIR="${BASE_DIR}/Validation_NCBI_10"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
OUTPUT_PARENT="${BASE_DIR}/nn_sensitivity_${RUN_TAG}"
mkdir -p "${OUTPUT_PARENT}"

NN_MODEL="${NN_MODEL:-${BASE_DIR}/slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet}"
NN_PREFILTER="${NN_PREFILTER:-0.0}"
NN_MULT_MODE="${NN_MULT_MODE:-ratio}"

# Default sweep ranges (override via env vars or NN_STRENGTH_VALUES / NN_CUTOFF_VALUES)
if [[ -n "${NN_STRENGTH_VALUES:-}" ]]; then
    read -ra STRENGTH_VALUES <<< "${NN_STRENGTH_VALUES}"
else
    STRENGTH_START="${NN_STRENGTH_START:-0}"
    STRENGTH_END="${NN_STRENGTH_END:-5}"
    STRENGTH_STEP="${NN_STRENGTH_STEP:-1}"
    mapfile -t STRENGTH_VALUES < <(seq "${STRENGTH_START}" "${STRENGTH_STEP}" "${STRENGTH_END}")
fi

if [[ -n "${NN_CUTOFF_VALUES:-}" ]]; then
    read -ra CUTOFF_VALUES <<< "${NN_CUTOFF_VALUES}"
else
    CUTOFF_START="${NN_CUTOFF_START:-0.3}"
    CUTOFF_END="${NN_CUTOFF_END:-0.9}"
    CUTOFF_STEP="${NN_CUTOFF_STEP:-0.1}"
    mapfile -t CUTOFF_VALUES < <(seq "${CUTOFF_START}" "${CUTOFF_STEP}" "${CUTOFF_END}")
fi

if [[ "${#STRENGTH_VALUES[@]}" -eq 0 || "${#CUTOFF_VALUES[@]}" -eq 0 ]]; then
    echo "No sweep values detected. Check NN_STRENGTH_* and NN_CUTOFF_* variables." >&2
    exit 1
fi

# Enforce CallGenes nnStrength constraint (0-1 range)
VALID_STRENGTH_VALUES=()
SKIPPED_STRENGTH_VALUES=()
for val in "${STRENGTH_VALUES[@]}"; do
    if (( $(echo "$val <= 1.0" | bc -l) )); then
        VALID_STRENGTH_VALUES+=("$val")
    else
        SKIPPED_STRENGTH_VALUES+=("$val")
    fi
done

if [[ "${#SKIPPED_STRENGTH_VALUES[@]}" -gt 0 ]]; then
    echo "[WARN] CallGenes enforces nnStrength <= 1. Skipping strengths: ${SKIPPED_STRENGTH_VALUES[*]}" >&2
fi

if [[ "${#VALID_STRENGTH_VALUES[@]}" -eq 0 ]]; then
    echo "[ERROR] All requested strengths exceed the allowed range (<=1). Aborting." >&2
    exit 1
fi

STRENGTH_VALUES=("${VALID_STRENGTH_VALUES[@]}")

shopt -s nullglob
GENOME_FILES=("${INPUT_DIR}"/*.fna.gz)
shopt -u nullglob

if [[ -n "${GENOME_LIMIT:-}" ]]; then
    if [[ "${GENOME_LIMIT}" =~ ^[0-9]+$ && "${GENOME_LIMIT}" -gt 0 ]]; then
        GENOME_FILES=("${GENOME_FILES[@]:0:${GENOME_LIMIT}}")
    else
        echo "[ERROR] GENOME_LIMIT must be a positive integer" >&2
        exit 1
    fi
fi

if [[ "${#GENOME_FILES[@]}" -eq 0 ]]; then
    echo "No .fna.gz inputs found in ${INPUT_DIR}" >&2
    exit 1
fi

PASS_MODES=("single_pass_nn" "two_pass_nn")
declare -A PASS_ARGS=(
    ["single_pass_nn"]=""
    ["two_pass_nn"]="2pass=t"
)

SUMMARY_CSV="${OUTPUT_PARENT}/summary_nn_sensitivity.csv"
echo "Genome,PassMode,NNStrength,NNCutoff,TP_Start,TP_Stop,FP_Start,FP_Stop,FN_Start,FN_Stop" > "${SUMMARY_CSV}"

TOTAL_RUNS=$(( ${#GENOME_FILES[@]} * ${#STRENGTH_VALUES[@]} * ${#CUTOFF_VALUES[@]} * ${#PASS_MODES[@]} ))
if [[ "${TOTAL_RUNS}" -eq 0 ]]; then
    echo "No runs scheduled. Check configured ranges and input set." >&2
    exit 1
fi

echo "Starting NN sensitivity sweep"
echo "  Genomes       : ${#GENOME_FILES[@]}"
echo "  NN strengths  : ${STRENGTH_VALUES[*]}"
echo "  NN cutoffs    : ${CUTOFF_VALUES[*]}"
echo "  Pass modes    : ${PASS_MODES[*]}"
echo "  Total run count: ${TOTAL_RUNS}"
echo "Using CURRENT_DIR=${CURRENT_DIR}"

RUN_INDEX=0

if ! pushd "${CURRENT_DIR}" >/dev/null; then
    echo "[ERROR] Unable to enter ${CURRENT_DIR}" >&2
    exit 1
fi

run_and_collect() {
    local genome_file="$1"
    local truth_gff="$2"
    local pass_mode="$3"
    local strength="$4"
    local cutoff="$5"

    local genome_name
    genome_name="$(basename "${genome_file}" .fna.gz)"
    local out_tag="${genome_name}.${pass_mode}.strength_${strength}_cutoff_${cutoff}"
    local out_gff="${OUTPUT_PARENT}/${out_tag}.gff"

    local args="${PASS_ARGS[$pass_mode]} nn=${NN_MODEL} strength=${strength} cutoff=${cutoff} prefilter_cutoff=${NN_PREFILTER} nnmultiplier_mode=${NN_MULT_MODE}"

    local callgenes_log="${OUTPUT_PARENT}/${out_tag}.callgenes.log"
    local compare_log="${OUTPUT_PARENT}/${out_tag}.compare.log"

    if ! java -ea -cp "." prok.CallGenes in="${genome_file}" outgff="${out_gff}" ${args} truegenes="${truth_gff}" log=t >"${callgenes_log}" 2>&1; then
        echo "[ERROR] CallGenes failed for ${genome_name} (${pass_mode}, strength=${strength}, cutoff=${cutoff}). See ${callgenes_log}" >&2
        return 1
    fi

    local metrics
    if ! java -ea -cp "." gff.CompareGff in="${out_gff}" ref="${truth_gff}" >"${compare_log}" 2>&1; then
        echo "[ERROR] CompareGff failed for ${genome_name} (${pass_mode}, strength=${strength}, cutoff=${cutoff}). See ${compare_log}" >&2
        return 1
    fi

    metrics=$(awk '
        /True Positive Start:/ {tp_s=$4}
        /True Positive Stop:/ {tp_st=$4}
        /False Positive Start:/ {fp_s=$4}
        /False Positive Stop:/ {fp_st=$4}
        /False Negative Start:/ {fn_s=$4}
        /False Negative Stop:/ {fn_st=$4}
        END {
            gsub(/,/, "", tp_s); gsub(/,/, "", tp_st);
            gsub(/,/, "", fp_s); gsub(/,/, "", fp_st);
            gsub(/,/, "", fn_s); gsub(/,/, "", fn_st);
            printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st
        }' "${compare_log}")

    echo "${genome_name},${pass_mode},${strength},${cutoff},${metrics}" >> "${SUMMARY_CSV}"
    rm -f "${out_gff}"
    rm -f "${callgenes_log}" "${compare_log}"
}

for cutoff in "${CUTOFF_VALUES[@]}"; do
    for strength in "${STRENGTH_VALUES[@]}"; do
        for genome_file in "${GENOME_FILES[@]}"; do
            truth_gff="${genome_file%.fna.gz}.gff.gz"
            if [[ ! -f "${truth_gff}" ]]; then
                echo "Skipping ${genome_file} (truth file ${truth_gff} missing)" >&2
                continue
            fi
            for pass_mode in "${PASS_MODES[@]}"; do
                ((RUN_INDEX+=1))
                echo "[${RUN_INDEX}/${TOTAL_RUNS}] Genome=${genome_file##*/} Pass=${pass_mode} Strength=${strength} Cutoff=${cutoff}"
                run_and_collect "${genome_file}" "${truth_gff}" "${pass_mode}" "${strength}" "${cutoff}"
            done
        done
    done
done

popd >/dev/null

if [[ "${RUN_INDEX}" -eq 0 ]]; then
    echo "[WARNING] No runs were executed. Check input/truth files and configuration." >&2
else
    echo "Sensitivity sweep complete."
    echo "Results stored in ${SUMMARY_CSV}"
fi
