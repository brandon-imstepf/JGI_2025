#!/bin/bash

usage(){
echo "
Written by Brandon Imstepf, with inspiration from Brian Bushnell.
Last modified June 30, 2025

Description: 
A main pipeline that processes an entire folder of genomes to create a
single, comprehensive machine learning training set. For each genome, it
generates true positives and high-confidence negatives, then concatenates
all results into a final TSV file.

Usage:   buildtrainingsetfolder.sh in=<folder> out=<master.tsv> [options]

Required Parameters:
in=<folder>     Input folder containing matching pairs of .fna.gz and .gff.gz files.
out=<master.tsv>  The final, combined output TSV training file.

Optional Pass-through Options (for CallGenes):
nofilter       Generates ALL candidates from CallGenes, not just filtered ones.
minlen=60      Minimum gene length for CallGenes predictions.
...any other valid CallGenes option...
"
}

# --- Boilerplate Setup ---
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null
CP="$DIR""current/"
z="-Xmx6g"
z2="-Xms6g"
# --- End Boilerplate ---

# --- Argument Parsing ---
INFOLDER=""
OUT_TSV=""
PASSTHRU_ARGS=()

for arg in "$@"; do
    case $arg in
        in=*) INFOLDER="${arg#in=}";;
        out=*) OUT_TSV="${arg#out=}";;
        -h|--help) usage; exit 0;;
        *) PASSTHRU_ARGS+=("$arg");;
    esac
done

if [[ -z "$INFOLDER" || -z "$OUT_TSV" ]]; then
    echo "ERROR: Missing required parameters: in=<folder> and out=<master.tsv>" >&2
    usage
    exit 1
fi

# --- Main Pipeline ---

TMPDIR=$(mktemp -d)
echo "Master pipeline started. Using temporary directory: $TMPDIR"


for ref_fasta in "$INFOLDER"/*.fna.gz; do
    base=$(basename "$ref_fasta" .fna.gz)
    ref_gff="$INFOLDER/$base.gff.gz"

    if [[ ! -f "$ref_gff" ]]; then
        echo "WARNING: No matching reference GFF for $ref_fasta. Skipping."
        continue
    fi

    echo "--- Processing Genome: $base ---"

    # Define temporary files for THIS specific genome
    ALL_GFF="$TMPDIR/${base}.all.gff"
    PRED_GFF="$TMPDIR/${base}.predicted.gff"
    UNION_GFF="$TMPDIR/${base}.union.gff"
    NEGATIVES_GFF="$TMPDIR/${base}.negatives.gff"
    POSITIVES_GFF="$TMPDIR/${base}.positives.gff"
    
    # STEP 1: Generate source GFFs
    echo "  [1/3] Generating gene calls..."
    # Get ALL candidates using 'nofilter'
    callgenes.sh in="$ref_fasta" outgff="$ALL_GFF" truegenes="$ref_gff" cds nofilter "${PASSTHRU_ARGS[@]}"
    # Get filtered, high-confidence predictions
    callgenes.sh in="$ref_fasta" outgff="$PRED_GFF" truegenes="$ref_gff" cds "${PASSTHRU_ARGS[@]}"

    # STEP 2: Perform set operations
    echo "  [2/3] Performing set operations..."
    # Union of predicted and reference
    gffsetop.sh in_a="$PRED_GFF" in_b="$ref_gff" out="$UNION_GFF" op=union
    # High-confidence negatives = all candidates minus the union
    gffsetop.sh in_a="$ALL_GFF" in_b="$UNION_GFF" out="$NEGATIVES_GFF" op=subtract
    # True positives = intersection of predicted and reference
    gffsetop.sh in_a="$PRED_GFF" in_b="$ref_gff" out="$POSITIVES_GFF" op=intersect

    # STEP 3: Convert to labeled TSV
    echo "  [3/3] Converting to labeled TSV format..."
    # Diagnostic: report GFF sizes and VECTOR-containing lines before conversion
    echo "    Negatives GFF lines: $(wc -l < "$NEGATIVES_GFF" 2>/dev/null || echo 0)"
    echo "    Negatives GFF VECTOR-containing lines: $(grep -c 'VECTOR=' "$NEGATIVES_GFF" 2>/dev/null || echo 0)"
    echo "    Positives GFF lines: $(wc -l < "$POSITIVES_GFF" 2>/dev/null || echo 0)"
    echo "    Positives GFF VECTOR-containing lines: $(grep -c 'VECTOR=' "$POSITIVES_GFF" 2>/dev/null || echo 0)"

    # Convert negatives with label 0, save to a final holding location
    gff2tsv.sh in="$NEGATIVES_GFF" out="$TMPDIR/${base}.negatives.tsv" label=0
    # Convert positives with label 1, save to a final holding location
    gff2tsv.sh in="$POSITIVES_GFF" out="$TMPDIR/${base}.positives.tsv" label=1

    # Diagnostic: report TSV sizes after conversion
    echo "    Negatives TSV lines: $(wc -l < "$TMPDIR/${base}.negatives.tsv" 2>/dev/null || echo 0)"
    echo "    Positives TSV lines: $(wc -l < "$TMPDIR/${base}.positives.tsv" 2>/dev/null || echo 0)"

    echo "--- Finished processing $base ---"
done

# --- Final Assembly ---
echo "Assembling final master training set: $OUT_TSV"

# Clear the output file and add a header
# Note: The header should match the columns from your GffToTsv
#echo -e "log_len\tstart_score\tmid_score\tend_score\tgene_gc\tgene_entropy\tcontig_gc\tcontig_entropy\t...\tmatch" > "$OUT_TSV"


# Balance operations
# Count total positives and negatives
POS_COUNT=$(find "$TMPDIR" -name '*.positives.tsv' -exec cat {} + | wc -l)
echo "Total positive examples: $POS_COUNT"
NEG_COUNT=$(find "$TMPDIR" -name '*.negatives.tsv' -exec cat {} + | wc -l)
echo "Total negative examples: $NEG_COUNT"

# Calculate max allowed negatives (2x positives)
MAX_NEG=$((POS_COUNT * 2))
echo "Maximum allowed negatives after balancing: $MAX_NEG"

# If too many negatives, subsample them
if (( NEG_COUNT > MAX_NEG )); then
    echo "Balancing negatives: subsampling $NEG_COUNT to $MAX_NEG"
    # First concatenate to a temp file
    find "$TMPDIR" -name '*.negatives.tsv' -exec cat {} + | shuf -n "$MAX_NEG" > "$TMPDIR/temp_balanced.tsv"
    # Then delete the originals
    find "$TMPDIR" -name '*.negatives.tsv' -delete
    # Finally rename to the proper pattern
    mv "$TMPDIR/temp_balanced.tsv" "$TMPDIR/balanced.negatives.tsv"
fi


# Concatenate all the individual positive and negative TSV files
# Use 'find' to gracefully handle cases where no files of a type were created
find "$TMPDIR" -name '*.positives.tsv' -print0 | xargs -0 cat >> "$OUT_TSV"
find "$TMPDIR" -name '*.negatives.tsv' -print0 | xargs -0 cat >> "$OUT_TSV"

# --- Cleanup ---
rm -r "$TMPDIR"
echo "Pipeline complete."
