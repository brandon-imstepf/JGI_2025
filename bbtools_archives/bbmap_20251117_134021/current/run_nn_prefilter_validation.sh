#!/bin/bash
set -euo pipefail

# This script benchmarks CallGenes NN runs with different prefilter_cutoff values.
# It iterates over Validation_NCBI_10, captures GFF + logs, and aggregates CompareGff metrics.

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INPUT_DIR="${BASE_DIR}/Validation_NCBI_10"
OUTPUT_DIR="${BASE_DIR}/nn_prefilter_validation"
NN_MODEL="${BASE_DIR}/slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_11-12_new.bbnet"
NNCUTOFF=0.6
NNSTRENGTH=1.0
PREFILTER_VALUES=("0.1" "0.0")
SUMMARY_FILE="${OUTPUT_DIR}/summary_prefilter.csv"
CALLGENES_SH="${BASE_DIR}/bbmap/callgenes.sh"
CLASSPATH="${BASE_DIR}/bbmap/current"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Validation folder not found at $INPUT_DIR" >&2
    exit 1
fi
if [[ ! -f "$NN_MODEL" ]]; then
    echo "ERROR: NN model not found at $NN_MODEL" >&2
    exit 1
fi
if [[ ! -x "$CALLGENES_SH" ]]; then
    echo "ERROR: callgenes.sh not executable at $CALLGENES_SH" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "Prefilter,Genome,TP_Start,TP_Stop,FP_Start,FP_Stop,FN_Start,FN_Stop" > "$SUMMARY_FILE"

for PREFILTER in "${PREFILTER_VALUES[@]}"; do
    RUN_DIR="${OUTPUT_DIR}/prefilter_${PREFILTER}"
    mkdir -p "$RUN_DIR"
    echo "=== Running validation with prefilter_cutoff=${PREFILTER} ==="

    for GENOME in "${INPUT_DIR}"/*.fna.gz; do
        BASE_NAME=$(basename "$GENOME" .fna.gz)
        TRUTH_GFF="${INPUT_DIR}/${BASE_NAME}.gff.gz"

        if [[ ! -f "$TRUTH_GFF" ]]; then
            echo "WARNING: Missing truth GFF for ${BASE_NAME}, skipping." >&2
            continue
        fi

        OUT_GFF="${RUN_DIR}/${BASE_NAME}.gff"
        STDOUT_LOG="${RUN_DIR}/${BASE_NAME}.stdout.log"

        echo "--- Processing ${BASE_NAME} (prefilter ${PREFILTER}) ---"
        "$CALLGENES_SH" \
            in="$GENOME" \
            outgff="$OUT_GFF" \
            compareto="$TRUTH_GFF" \
            truegenes="$TRUTH_GFF" \
            nn="$NN_MODEL" \
            nncutoff="$NNCUTOFF" \
            strength="$NNSTRENGTH" \
            prefilter_cutoff="$PREFILTER" \
            log=t > "$STDOUT_LOG" 2>&1

        METRICS=$(java -ea -cp "$CLASSPATH" gff.CompareGff in="$OUT_GFF" ref="$TRUTH_GFF" 2>&1 | awk '
            /True Positive Start:/ {tp_s=$4}
            /True Positive Stop:/  {tp_st=$4}
            /False Positive Start:/ {fp_s=$4}
            /False Positive Stop:/ {fp_st=$4}
            /False Negative Start:/ {fn_s=$4}
            /False Negative Stop:/ {fn_st=$4}
            END {
                gsub(/,/, "", tp_s); gsub(/,/, "", tp_st);
                gsub(/,/, "", fp_s); gsub(/,/, "", fp_st);
                gsub(/,/, "", fn_s); gsub(/,/, "", fn_st);
                printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st
            }')

        echo "${PREFILTER},${BASE_NAME},${METRICS}" >> "$SUMMARY_FILE"
    done
done

echo "Prefilter benchmarking complete. Summary saved to ${SUMMARY_FILE}"
