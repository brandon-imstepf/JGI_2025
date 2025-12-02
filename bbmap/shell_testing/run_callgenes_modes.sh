#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BBMAP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CURRENT_DIR="${BBMAP_DIR}/current"
BASE_DIR="$(cd "${BBMAP_DIR}/.." && pwd)"
INPUT_DIR="${BASE_DIR}/Validation_NCBI_10"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
OUTPUT_PARENT="${BASE_DIR}/callgenes_modes_${RUN_TAG}"
NN_MODEL="${BASE_DIR}/slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet" # Path to the neural network model
NN_STRENGTH=1.0
NN_CUTOFF=0.6
NN_PREFILTER=0.0
NN_MULT_MODE="${NN_MULT_MODE:-ratio}"
NN_MODE_LABEL="${NN_MULT_MODE}" # may be empty when no NN
if [[ -z "$NN_MODE_LABEL" || "$NN_MODE_LABEL" == "ratio" && -z "${NN_MULT_MODE+x}" ]]; then
    NN_MODE_LABEL="baseline"
fi
OUTPUT_DIR="${OUTPUT_PARENT}/${NN_MODE_LABEL}"

# --- Setup ---
mkdir -p "$OUTPUT_DIR"
echo "Mode,Genome,TP_Start,TP_Stop,FP_Start,FP_Stop,FN_Start,FN_Stop" > "${OUTPUT_DIR}/summary_results_callgenes.csv"

pushd "$CURRENT_DIR" >/dev/null

# --- Main Loop ---
for GENOME_FILE in "${INPUT_DIR}"/*.fna.gz; do
    BASE_NAME=$(basename "$GENOME_FILE" .fna.gz)
    TRUTH_GFF="${INPUT_DIR}/${BASE_NAME}.gff.gz"

    echo "--- Processing: $BASE_NAME ---"

    # --- Define a function for running and evaluating a mode ---
    run_and_evaluate() {
        local MODE_NAME="$1"
        local CALLGENES_ARGS="$2"
        
        echo "  Running Mode: $MODE_NAME..."
        local OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.${MODE_NAME}.gff"
        local CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.${MODE_NAME}.cmp"

        # Run CallGenes with the specified arguments
        java -ea -cp "." prok.CallGenes in="$GENOME_FILE" outgff="$OUT_GFF" ${CALLGENES_ARGS} truegenes="$TRUTH_GFF" log=t > /dev/null 2>&1

        # Run CompareGff and parse metrics in one step
        local METRICS=$(java -ea -cp "." gff.CompareGff in="$OUT_GFF" ref="$TRUTH_GFF" 2>&1 | awk ' \
            /True Positive Start:/ {tp_s=$4} \
            /True Positive Stop:/ {tp_st=$4} \
            /False Positive Start:/ {fp_s=$4} \
            /False Positive Stop:/ {fp_st=$4} \
            /False Negative Start:/ {fn_s=$4} \
            /False Negative Stop:/ {fn_st=$4} \
            END {gsub(/,/, "", tp_s); gsub(/,/, "", tp_st); gsub(/,/, "", fp_s); gsub(/,/, "", fp_st); gsub(/,/, "", fn_s); gsub(/,/, "", fn_st); printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st}')
        
        # Append to summary CSV
        echo "$MODE_NAME,$BASE_NAME,$METRICS" >> "${OUTPUT_DIR}/summary_results_callgenes.csv"
        
        # Cleanup intermediate files for this run
        rm "$OUT_GFF"
    }

    # --- Run all four modes ---
    run_and_evaluate "single_pass" ""
    run_and_evaluate "two_pass" "2pass=t"
    run_and_evaluate "single_pass_nn" "nn=${NN_MODEL} strength=${NN_STRENGTH} cutoff=${NN_CUTOFF} prefilter_cutoff=${NN_PREFILTER} nnmultiplier_mode=${NN_MULT_MODE}"
    run_and_evaluate "two_pass_nn" "2pass=t nn=${NN_MODEL} strength=${NN_STRENGTH} cutoff=${NN_CUTOFF} prefilter_cutoff=${NN_PREFILTER} nnmultiplier_mode=${NN_MULT_MODE}"

done

echo "--- All CallGenes Mode Benchmarking Complete ---"
echo "Results aggregated in ${OUTPUT_DIR}/summary_results_callgenes.csv"

popd >/dev/null
