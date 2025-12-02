#!/bin/bash
set -e

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BBMAP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CURRENT_DIR="${BBMAP_DIR}/current"
BASE_DIR="$(cd "${BBMAP_DIR}/.." && pwd)"
INPUT_DIR="${BASE_DIR}/Validation_NCBI_10"
OUTPUT_DIR="${BASE_DIR}/benchmarking_results_all"
PRODIGAL_EXEC="${BASE_DIR}/comparison_models/Prodigal/prodigal"
GLIMMER_DIR="${BASE_DIR}/comparison_models/glimmer/glimmer3.02"
GENEMARK_EXEC="${BASE_DIR}/comparison_models/genemark/gms2_linux_64/gms2.pl"
CALLGENES_EXEC="${BBMAP_DIR}/callgenes.sh"
NN_FILE="${BASE_DIR}/slurm/cycles_trained_8dims_400k_ncbi_total_8f_6l_small_11-12_new.bbnet"
NN_STRENGTH=0.5
NN_CUTOFF=0.5
NN_PREFILTER=0.1

collect_metrics() {
    local tool_label="$1"
    local prediction_gff="$2"
    local truth_gff="$3"
    local genome_name="$4"

    local metrics
    metrics=$(java -ea -cp "." gff.CompareGff in="$prediction_gff" ref="$truth_gff" 2>&1 | awk ' \
        /True Positive Start:/ {tp_s=$4} \
        /True Positive Stop:/ {tp_st=$4} \
        /False Positive Start:/ {fp_s=$4} \
        /False Positive Stop:/ {fp_st=$4} \
        /False Negative Start:/ {fn_s=$4} \
        /False Negative Stop:/ {fn_st=$4} \
        END {gsub(/,/, "", tp_s); gsub(/,/, "", tp_st); gsub(/,/, "", fp_s); gsub(/,/, "", fp_st); gsub(/,/, "", fn_s); gsub(/,/, "", fn_st); printf "%s,%s,%s,%s,%s,%s", tp_s, tp_st, fp_s, fp_st, fn_s, fn_st}')

    echo "${tool_label},${genome_name},${metrics}" >> "${OUTPUT_DIR}/summary_results.csv"
}

# --- Setup ---
mkdir -p "$OUTPUT_DIR"
echo "Tool,Genome,TP_Start,TP_Stop,FP_Start,FP_Stop,FN_Start,FN_Stop" > "${OUTPUT_DIR}/summary_results.csv"

pushd "$CURRENT_DIR" >/dev/null

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
    collect_metrics "Prodigal" "$PRODIGAL_OUT_GFF" "$TRUTH_GFF" "$BASE_NAME"

    # --- 2. Run Glimmer ---
    echo "  Running Glimmer..."
    GLIMMER_SCRIPT_MODIFIED="${OUTPUT_DIR}/g3-iterated-modified.csh"
    GLIMMER_TAG="${OUTPUT_DIR}/${BASE_NAME}"
    GLIMMER_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.glimmer.gff"

    # Create (or refresh) a modified Glimmer driver that uses the bundled binaries
    cat << EOF > "$GLIMMER_SCRIPT_MODIFIED"
#!/bin/csh
set genome = \$1
set tag = \$2
set awkpath = "${GLIMMER_DIR}/scripts"
set glimmerpath = "${GLIMMER_DIR}/bin"
set glimmeropts = "-o50 -g110 -t30"
\$glimmerpath/long-orfs -n -t 1.15 \$genome \$tag.longorfs >& /dev/null
\$glimmerpath/extract -t \$genome \$tag.longorfs > \$tag.train
\$glimmerpath/build-icm -r \$tag.icm < \$tag.train >& /dev/null
\$glimmerpath/glimmer3 \$glimmeropts \$genome \$tag.icm \$tag.run1 >& /dev/null
tail -n +2 \$tag.run1.predict > \$tag.coords
set startuse = \`\$glimmerpath/start-codon-distrib -3 \$genome \$tag.coords\`
\$glimmerpath/glimmer3 \$glimmeropts -P \$startuse \$genome \$tag.icm \$tag >& /dev/null
EOF
    chmod +x "$GLIMMER_SCRIPT_MODIFIED"

    (cd "$OUTPUT_DIR" && ./"$(basename "$GLIMMER_SCRIPT_MODIFIED")" "${BASE_NAME}.fna" "$BASE_NAME") > /dev/null 2>&1

    # Convert the predict file to GFF, normalizing coordinate order for reverse-strand calls
    grep -v '^>' "${GLIMMER_TAG}.predict" | grep -v '^$' | \
        awk -v SEQ_ID="$SEQ_ID" 'BEGIN {OFS="\t"} {
            start=$2; end=$3;
            if (start > end) { tmp=start; start=end; end=tmp; }
            strand = ($4 ~ /^-/ ? "-" : "+");
            print SEQ_ID, "Glimmer3", "CDS", start, end, $5, strand, "0", "gene_id \""$1"\";";
        }' > "$GLIMMER_OUT_GFF"

    collect_metrics "Glimmer" "$GLIMMER_OUT_GFF" "$TRUTH_GFF" "$BASE_NAME"

    # --- 3. Run GeneMarkS2 ---
    echo "  Running GeneMarkS2..."
    GENEMARK_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.genemark.gff"
    GENEMARK_CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.genemark.cmp"

    perl "$GENEMARK_EXEC" --seq "${OUTPUT_DIR}/${BASE_NAME}.fna" --genome-type bacteria --format gff3 --output "$GENEMARK_OUT_GFF"
    GENEMARK_CDS_GFF="${GENEMARK_OUT_GFF%.gff}.cds.gff"
    awk '$3=="CDS"' "$GENEMARK_OUT_GFF" > "$GENEMARK_CDS_GFF"
    collect_metrics "GeneMarkS2" "$GENEMARK_CDS_GFF" "$TRUTH_GFF" "$BASE_NAME"

    # --- 4. Run CallGenes baseline (1-pass) ---
    echo "  Running CallGenes baseline (1-pass)..."
    CALLGENES_SINGLE_OUT="${OUTPUT_DIR}/${BASE_NAME}.cg.single.gff"
    "$CALLGENES_EXEC" in="${OUTPUT_DIR}/${BASE_NAME}.fna" outgff="$CALLGENES_SINGLE_OUT" passes=1 > /dev/null 2>&1
    collect_metrics "CallGenes-1pass" "$CALLGENES_SINGLE_OUT" "$TRUTH_GFF" "$BASE_NAME"

    # --- 5. Run CallGenes baseline (2-pass) ---
    echo "  Running CallGenes baseline (2-pass)..."
    CALLGENES_DOUBLE_OUT="${OUTPUT_DIR}/${BASE_NAME}.cg.double.gff"
    "$CALLGENES_EXEC" in="${OUTPUT_DIR}/${BASE_NAME}.fna" outgff="$CALLGENES_DOUBLE_OUT" passes=2 > /dev/null 2>&1
    collect_metrics "CallGenes-2pass" "$CALLGENES_DOUBLE_OUT" "$TRUTH_GFF" "$BASE_NAME"

    # --- 6. Run CallGenes NN (1-pass) ---
    echo "  Running CallGenes NN (1-pass)..."
    CALLGENES_NN_SINGLE_OUT="${OUTPUT_DIR}/${BASE_NAME}.cg.nn.single.gff"
    "$CALLGENES_EXEC" in="${OUTPUT_DIR}/${BASE_NAME}.fna" outgff="$CALLGENES_NN_SINGLE_OUT" passes=1 net="$NN_FILE" strength="$NN_STRENGTH" cutoff="$NN_CUTOFF" prefilter_cutoff="$NN_PREFILTER" > /dev/null 2>&1
    collect_metrics "CallGenes-NN-1pass" "$CALLGENES_NN_SINGLE_OUT" "$TRUTH_GFF" "$BASE_NAME"

    # --- 7. Run CallGenes NN (2-pass) ---
    echo "  Running CallGenes NN (2-pass)..."
    CALLGENES_NN_DOUBLE_OUT="${OUTPUT_DIR}/${BASE_NAME}.cg.nn.double.gff"
    "$CALLGENES_EXEC" in="${OUTPUT_DIR}/${BASE_NAME}.fna" outgff="$CALLGENES_NN_DOUBLE_OUT" passes=2 net="$NN_FILE" strength="$NN_STRENGTH" cutoff="$NN_CUTOFF" prefilter_cutoff="$NN_PREFILTER" > /dev/null 2>&1
    collect_metrics "CallGenes-NN-2pass" "$CALLGENES_NN_DOUBLE_OUT" "$TRUTH_GFF" "$BASE_NAME"

    # --- Cleanup ---
    rm "${OUTPUT_DIR}/${BASE_NAME}.fna"
done

echo "--- All Benchmarking Complete ---"
echo "Results aggregated in ${OUTPUT_DIR}/summary_results.csv"

popd >/dev/null
