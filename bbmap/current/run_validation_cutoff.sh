#!/bin/bash
set -e # Exit on any error

# Change to the script's directory to ensure relative paths work correctly
cd "$(dirname "$0")/../.."

# --- Configuration ---
# This script should be run from the project root directory that contains 'bbmap', 'Validation_NCBI_10', etc.
INPUT_FOLDER="./Validation_NCBI_10"
BASE_OUTPUT_DIR="./nn_cutoff_validation"
NNCUTOFFS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
PREFILTER_CUTOFF=0.15
STRENGTH=1.0
NN_PATH="./slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet"
RESULTS_FILE="${BASE_OUTPUT_DIR}/validation_results.csv"

# --- Setup ---
mkdir -p "$BASE_OUTPUT_DIR"
echo "NNCutoff,TotalExecutionTime,TotalTruePositives,TotalFalsePositives,TotalFalseNegatives" > "$RESULTS_FILE"

# --- Outer Loop: Iterate over each nncutoff value ---
for nncutoff in "${NNCUTOFFS[@]}"; do
    echo "=================================================="
    echo "RUNNING VALIDATION WITH NN CUTOFF: $nncutoff"
    echo "=================================================="

    # Directory to store all outputs for this specific nncutoff run
    NNCUTOFF_OUTPUT_DIR="${BASE_OUTPUT_DIR}/nncutoff_${nncutoff}"
    mkdir -p "$NNCUTOFF_OUTPUT_DIR"

    # --- Inner Loop: Iterate over each genome in the input folder ---
    for fna_file in "${INPUT_FOLDER}"/*.fna.gz; do
        
        base=$(basename "$fna_file" .fna.gz)
        echo "--- Processing: $base ---"

        gff_file="${INPUT_FOLDER}/${base}.gff.gz"
        if [ ! -f "$gff_file" ]; then
            echo "WARNING: Corresponding GFF file not found for $fna_file. Skipping."
            continue
        fi

        output_gff="${NNCUTOFF_OUTPUT_DIR}/${base}.gff"
        time_log="${NNCUTOFF_OUTPUT_DIR}/${base}.time.log"

        # Run the standard callgenes.sh script for a single genome
        {
            time -p ./bbmap/callgenes.sh \
                in="$fna_file" \
                outgff="$output_gff" \
                compareto="$gff_file" \
                truegenes="$gff_file" \
                prefilter_cutoff="$PREFILTER_CUTOFF" \
                strength="$STRENGTH" \
                nncutoff="$nncutoff" \
                nn="$NN_PATH" \
                log=t
        } 2> "$time_log"
    done

    # --- Aggregation Step for the current nncutoff ---
    echo "--- Aggregating results for nncutoff $nncutoff ---"
    
    TOTAL_TP=0
    TOTAL_FP=0
    TOTAL_FN=0
    TOTAL_TIME=0.0

    for base_gff_file in "${NNCUTOFF_OUTPUT_DIR}"/*.gff; do
        base_name=$(basename "$base_gff_file" .gff)
        
        log_dir="${NNCUTOFF_OUTPUT_DIR}/${base_name}"
        log_file="${log_dir}/log.txt"
        time_file="${NNCUTOFF_OUTPUT_DIR}/${base_name}.time.log"

        if [ -f "$log_file" ]; then
            TP=$(grep "True positives:" "$log_file" | awk '{print $3}')
            FP=$(grep "False positives:" "$log_file" | awk '{print $3}')
            FN=$(grep "False negatives:" "$log_file" | awk '{print $3}')
            
            TOTAL_TP=$((TOTAL_TP + ${TP:-0}))
            TOTAL_FP=$((TOTAL_FP + ${FP:-0}))
            TOTAL_FN=$((TOTAL_FN + ${FN:-0}))
        else
            echo "WARNING: Log file not found for ${base_name}"
        fi

        if [ -f "$time_file" ]; then
            EXEC_TIME=$(grep "real" "$time_file" | awk '{print $2}')
            TOTAL_TIME=$(echo "$TOTAL_TIME + $EXEC_TIME" | bc)
        else
            echo "WARNING: Time log not found for ${base_name}"
        fi
    done

    # --- Record the final aggregated results for this nncutoff ---
    echo "Total Execution Time: ${TOTAL_TIME}s"
    echo "Total True Positives: $TOTAL_TP"
    echo "Total False Positives: $TOTAL_FP"
    echo "Total False Negatives: $TOTAL_FN"
    echo "${nncutoff},${TOTAL_TIME},${TOTAL_TP},${TOTAL_FP},${TOTAL_FN}" >> "$RESULTS_FILE"
    echo ""
done

echo "=================================================="
echo "All validation experiments complete."
echo "Aggregated results saved to ${RESULTS_FILE}"
echo "=================================================="