#!/bin/bash
set -e

# --- Configuration ---
BASE_DIR=../..
INPUT_DIR="${BASE_DIR}/sh_test"
OUTPUT_DIR="${BASE_DIR}/benchmarking_results"
PRODIGAL_EXEC="${BASE_DIR}/comparison_models/Prodigal/prodigal"
GLIMMER_DIR="$(pwd)/${BASE_DIR}/comparison_models/glimmer/glimmer3.02"

# --- Setup ---
mkdir -p "$OUTPUT_DIR"

# --- Run on a single file for testing ---
GENOME_FILE=$(find "$INPUT_DIR" -name "*.fna.gz" | head -n 1)
BASE_NAME=$(basename "$GENOME_FILE" .fna.gz)
TRUTH_GFF="${INPUT_DIR}/${BASE_NAME}.gff.gz"
SEQ_ID=$(gunzip -c "$GENOME_FILE" | head -n 1 | cut -d' ' -f1 | sed 's/>//')

echo "--- Processing: $BASE_NAME (SEQ_ID: $SEQ_ID) ---"

# --- 1. Run Prodigal ---
echo "Running Prodigal..."
PRODIGAL_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.prodigal.gff"
gunzip -c "$GENOME_FILE" > "${OUTPUT_DIR}/${BASE_NAME}.fna" # Prodigal doesn't take gzipped input

"$PRODIGAL_EXEC" -i "${OUTPUT_DIR}/${BASE_NAME}.fna" -o "$PRODIGAL_OUT_GFF" -f gff

echo "Prodigal complete. Output at: $PRODIGAL_OUT_GFF"
rm "${OUTPUT_DIR}/${BASE_NAME}.fna" # Clean up unzipped file

# --- 2. Run Glimmer ---
echo "Running Glimmer..."
GLIMMER_SCRIPT_MODIFIED="${OUTPUT_DIR}/g3-iterated-modified.csh"
GLIMMER_TAG="${OUTPUT_DIR}/${BASE_NAME}" # Used by Glimmer to prefix its output files

# Create a new, corrected glimmer script using a heredoc
cat << EOF > "$GLIMMER_SCRIPT_MODIFIED"
#!/bin/csh

if (\$#argv < 2) then
  echo "Usage:  g3-iterated.csh  <genome> <tag>"
  exit -1;
endif

set genome = \$1
set tag = \$2

# --- Corrected Paths ---
set awkpath = "${GLIMMER_DIR}/scripts"
set glimmerpath = "${GLIMMER_DIR}/bin"

# add/change glimmer options here
set glimmeropts = "-o50 -g110 -t30"

# Step 1: Find long, non-overlapping orfs to use as a training set
echo "Step 1 of 8:  Finding long orfs for training"
\$glimmerpath/long-orfs -n -t 1.15 \$genome \$tag.longorfs

# Step 2: Extract the training sequences from the genome file
echo "Step 2 of 8:  Extracting training sequences"
\$glimmerpath/extract -t \$genome \$tag.longorfs > \$tag.train

# Step 3: Build the icm from the training sequences
echo "Step 3 of 8:  Building ICM"
\$glimmerpath/build-icm -r \$tag.icm < \$tag.train

# Step 4: Run first Glimmer
echo "Step 4 of 8:  Running first Glimmer3"
\$glimmerpath/glimmer3 \$glimmeropts \$genome \$tag.icm \$tag.run1

# Step 5: Get training coordinates from first predictions
echo "Step 5 of 8:  Getting training coordinates"
tail -n +2 \$tag.run1.predict > \$tag.coords

# Step 6: Create a position weight matrix (PWM) - SKIPPED
echo "Step 6 of 8:  Making PWM from upstream regions (SKIPPED - ELPH not installed)"

# Step 7: Determine the distribution of start-codon usage
echo "Step 7 of 8:  Getting start-codon usage"
set startuse = \`\$glimmerpath/start-codon-distrib -3 \$genome \$tag.coords\`

# Step 8: Run second Glimmer
echo "Step 8 of 8:  Running second Glimmer3"
# Note: The -b (motif) option is removed because ELPH was not run
\$glimmerpath/glimmer3 \$glimmeropts -P \$startuse \$genome \$tag.icm \$tag

EOF

chmod +x "$GLIMMER_SCRIPT_MODIFIED"

# Glimmer also needs an unzipped genome
gunzip -c "$GENOME_FILE" > "${OUTPUT_DIR}/${BASE_NAME}.fna"

# Run the modified glimmer script
# The script takes the genome and an output tag as arguments
(cd "$OUTPUT_DIR" && ./"$(basename "$GLIMMER_SCRIPT_MODIFIED")" "${BASE_NAME}.fna" "$BASE_NAME")

GLIMMER_PREDICT_FILE="${GLIMMER_TAG}.predict"
echo "Glimmer prediction output at: $GLIMMER_PREDICT_FILE"

# Convert Glimmer's .predict file to GFF format.
GLIMMER_OUT_GFF="${OUTPUT_DIR}/${BASE_NAME}.glimmer.gff"
echo "Converting Glimmer output to GFF..."
# The awk script reformats the .predict file into GFF.
# It handles the strand format (+1 -> +) and creates a valid GFF structure.
awk -v SEQ_ID="$SEQ_ID" 'BEGIN {OFS="\t"} !/^>/ { \
    strand = ($4 ~ /^-/ ? "-" : "+"); \
    print SEQ_ID, "Glimmer", "CDS", $2, $3, $5, strand, ".", "gene_id \""$1"\";" \
}' "$GLIMMER_PREDICT_FILE" > "$GLIMMER_OUT_GFF"

echo "Glimmer GFF output at: $GLIMMER_OUT_GFF"

rm "${OUTPUT_DIR}/${BASE_NAME}.fna" # Clean up unzipped file

# --- 3. Compare Results ---
echo "Comparing results to truth GFF..."
PRODIGAL_CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.prodigal.cmp"
GLIMMER_CMP_OUT="${OUTPUT_DIR}/${BASE_NAME}.glimmer.cmp"

java -ea -cp "." gff.CompareGff in="$PRODIGAL_OUT_GFF" ref="$TRUTH_GFF" > "$PRODIGAL_CMP_OUT" 2>&1
java -ea -cp "." gff.CompareGff in="$GLIMMER_OUT_GFF" ref="$TRUTH_GFF" > "$GLIMMER_CMP_OUT" 2>&1

echo "Comparison complete. Results in:"
echo "  - $PRODIGAL_CMP_OUT"
echo "  - $GLIMMER_CMP_OUT"
