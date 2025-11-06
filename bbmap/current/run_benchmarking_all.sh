#!/bin/bash
set -e

# --- Configuration ---
BASE_DIR=../..
INPUT_DIR="${BASE_DIR}/Validation_NCBI_10"
OUTPUT_DIR="${BASE_DIR}/benchmarking_results_all"
PRODIGAL_EXEC="${BASE_DIR}/comparison_models/Prodigal/prodigal"
GLIMMER_DIR="$(pwd)/${BASE_DIR}/comparison_models/glimmer/glimmer3.02"
GENEMARK_EXEC="${BASE_DIR}/comparison_models/genemark/gms2_linux_64/gms2.pl"

# --- Setup ---
mkdir -p "$OUTPUT_DIR"
echo "Tool,Genome,TP_Start,TP_Stop,FP_Start,FP_Stop,FN_Start,FN_Stop" > "${OUTPUT_DIR}/summary_results.csv"

# --- Main Loop ---
for GENOME_FILE in "${INPUT_DIR}"/*.fna.gz; do
    BASE_NAME=$(basename "$GENOME_FILE" .fna.gz)
    TRUTH_GFF="${INPUT_DIR}/${BASE_NAME}.gff.gz"
    SEQ_ID=$(gunzip -c "$GENOME_FILE" | head -n 1 | cut -d' ' -f1 | sed 's/>//')

    echo "--- Processing: $BASE_NAME (SEQ_ID: $SEQ_ID) ---"

    # --- 1. Run Prodigal ---
    echo "  Running Prodigal..."
    PRODIGAL_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.prodigal.gff"
    gunzip -c "$GENOME_FILE" > "${OUTPUT_DIR}/${BASE_NAME}.fna"

    "$PRODIGAL_EXEC" -i "${OUTPUT_DIR}/${BASE_NAME}.fna" -o "$PRODIGAL_OUT_GFF" -f gff > /dev/null 2>&1
    
    PRODIGAL_METRICS=$(java -ea -cp "." gff.CompareGff in="$PRODIGAL_OUT_GFF" ref="$TRUTH_GFF" 2>&1 | awk ' \
        /True Positive Start:/ {tp_s=$4} \
        /True Positive Stop:/ {tp_st=$4} \
        /False Positive Start:/ {fp_s=$4} \
        /False Positive Stop:/ {fp_st=$4} \
        /False Negative Start:/ {fn_s=$4} \
        /False Negative Stop:/ {fn_st=$4} \
        END {gsub(/,/, "", tp_s); gsub(/,/, "", tp_st); gsub(/,/, "", fp_s); gsub(/,/, "", fp_st); gsub(/,/, "", fn_s); gsub(/,/, "", fn_st); printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st}')
    echo "Prodigal,$BASE_NAME,$PRODIGAL_METRICS" >> "${OUTPUT_DIR}/summary_results.csv"

    # --- 2. Run Glimmer (Commented out) ---
    # echo "  Running Glimmer..."
    # GLIMMER_SCRIPT_MODIFIED="${OUTPUT_DIR}/g3-iterated-modified.csh"
    # GLIMMER_TAG="${OUTPUT_DIR}/${BASE_NAME}"
    # GLIMMER_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.glimmer.gff"
    #
    # # Create the modified glimmer script
    # # (This is repeated each loop, but is fast and ensures correctness)
    # cat << EOF > "$GLIMMER_SCRIPT_MODIFIED"
    # #!/bin/csh
    # set genome = \$1
    # set tag = \$2
    # set awkpath = "${GLIMMER_DIR}/scripts"
    # set glimmerpath = "${GLIMMER_DIR}/bin"
    # set glimmeropts = "-o50 -g110 -t30"
    # \$glimmerpath/long-orfs -n -t 1.15 \$genome \$tag.longorfs >& /dev/null
    # \$glimmerpath/extract -t \$genome \$tag.longorfs > \$tag.train
    # \$glimmerpath/build-icm -r \$tag.icm < \$tag.train >& /dev/null
    # \$glimmerpath/glimmer3 \$glimmeropts \$genome \$tag.icm \$tag.run1 >& /dev/null
    # tail -n +2 \$tag.run1.predict > \$tag.coords
    # set startuse = \`\$glimmerpath/start-codon-distrib -3 \$genome \$tag.coords\`
    # \$glimmerpath/glimmer3 \$glimmeropts -P \$startuse \$genome \$tag.icm \$tag >& /dev/null
    # EOF
    # chmod +x "$GLIMMER_SCRIPT_MODIFIED"
    #
    # (cd "$OUTPUT_DIR" && ./"$(basename "$GLIMMER_SCRIPT_MODIFIED")" "${BASE_NAME}.fna" "$BASE_NAME") > /dev/null 2>&1
    #
    # # The g3-iterated.csh script uses tail +2, but the .predict file often has a header.
    # # A more robust solution is to filter out any line that starts with '>' or is empty.
    # grep -v '^>' "${GLIMMER_TAG}.predict" | grep -v '^$' | awk -v SEQ_ID="$SEQ_ID" 'BEGIN {OFS="\t"} { \
    #     strand = ($4 ~ /^-/ ? "-" : "+"); \
    #     print SEQ_ID, "Glimmer3", "CDS", $2, $3, $5, strand, "0", "gene_id \""$1"\";"; \
    # }' > "$GLIMMER_OUT_GFF"
    #
    # GLIMMER_CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.glimmer.cmp"
    # java -ea -cp "." gff.CompareGff in="$GLIMMER_OUT_GFF" ref="$TRUTH_GFF" > "$GLIMMER_CMP_OUT" 2>&1
    #
    # TP_START=$(grep "True Positive Start:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # TP_STOP=$(grep "True Positive Stop:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # FP_START=$(grep "False Positive Start:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # FP_STOP=$(grep "False Positive Stop:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # FN_START=$(grep "False Negative Start:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # FN_STOP=$(grep "False Negative Stop:" "$GLIMMER_CMP_OUT" | awk '{print $4}')
    # echo "Glimmer,$BASE_NAME,$TP_START,$TP_STOP,$FP_START,$FP_STOP,$FN_START,$FN_STOP" >> "${OUTPUT_DIR}/summary_results.csv"

    # --- 3. Run GeneMarkS2 ---
    echo "  Running GeneMarkS2..."
    GENEMARK_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.genemark.gff"
    GENEMARK_CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.genemark.cmp"

    perl "$GENEMARK_EXEC" --seq "${OUTPUT_DIR}/${BASE_NAME}.fna" --genome-type bacteria --format gff3 --output "$GENEMARK_OUT_GFF"
    GENEMARK_METRICS=$(java -ea -cp "." gff.CompareGff in="$GENEMARK_OUT_GFF" ref="$TRUTH_GFF" 2>&1 | awk ' \
        /True Positive Start:/ {tp_s=$4} \
        /True Positive Stop:/ {tp_st=$4} \
        /False Positive Start:/ {fp_s=$4} \
        /False Positive Stop:/ {fp_st=$4} \
        /False Negative Start:/ {fn_s=$4} \
        /False Negative Stop:/ {fn_st=$4} \
        END {gsub(/,/, "", tp_s); gsub(/,/, "", tp_st); gsub(/,/, "", fp_s); gsub(/,/, "", fp_st); gsub(/,/, "", fn_s); gsub(/,/, "", fn_st); printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st}')
    echo "GeneMarkS2,$BASE_NAME,$GENEMARK_METRICS" >> "${OUTPUT_DIR}/summary_results.csv"

    # --- Cleanup ---
    rm "${OUTPUT_DIR}/${BASE_NAME}.fna"
done

echo "--- All Benchmarking Complete ---"
echo "Results aggregated in ${OUTPUT_DIR}/summary_results.csv"