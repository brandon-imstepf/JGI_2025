#!/bin/bash

usage(){
echo "
Description: 
Builds a complete ML training set from a reference genome and GFF.
This script generates true positives and high-confidence negatives.

Usage:   build_training_set.sh ref_fasta=<file> ref_gff=<file> out=<file.tsv>

Required Parameters:
ref_fasta=<file>  The reference genome FASTA file.
ref_gff=<file>    The reference (ground truth) GFF file.
out=<file.tsv>    The final, combined output TSV training file.
"
}

# --- Argument Parsing ---
REF_FASTA=""
REF_GFF=""
OUT_TSV=""

for arg in "$@"; do
    case $arg in
        ref_fasta=*) REF_FASTA="${arg#ref_fasta=}";;
        ref_gff=*) REF_GFF="${arg#ref_gff=}";;
        out=*) OUT_TSV="${arg#out=}";;
        -h|--help) usage; exit 0;;
    esac
done

if [[ -z "$REF_FASTA" || -z "$REF_GFF" || -z "$OUT_TSV" ]]; then
    echo "ERROR: All parameters are required." >&2
    usage
    exit 1
fi

# --- Main Pipeline ---
TMPDIR=$(mktemp -d)
echo "Pipeline started. Using temporary directory: $TMPDIR"

# --- 3.1: Generate source GFFs ---
echo "[1/4] Generating raw gene calls from FASTA..."
ALL_GFF="$TMPDIR/all_candidates.gff"
PRED_GFF="$TMPDIR/predicted.gff"

# 3.1a: Get ALL candidates
callgenes.sh in="$REF_FASTA" outgff="$ALL_GFF" nofilter
# 3.1b: Get filtered, high-confidence predictions
callgenes.sh in="$REF_FASTA" outgff="$PRED_GFF"

# --- 3.2: Perform set operations ---
echo "[2/4] Performing set operations to isolate feature sets..."
UNION_GFF="$TMPDIR/union.gff"
NEGATIVES_GFF="$TMPDIR/negatives.gff"
POSITIVES_GFF="$TMPDIR/positives.gff"

# 3.2a: Union of predicted and reference
gffsetop.sh in_a="$PRED_GFF" in_b="$REF_GFF" out="$UNION_GFF" op=union
# 3.2b: High-confidence negatives = all candidates minus the union
gffsetop.sh in_a="$ALL_GFF" in_b="$UNION_GFF" out="$NEGATIVES_GFF" op=subtract
# 3.2c: High confidence positives = intersection of predicted and reference
gffsetop.sh in_a="$PRED_GFF" in_b="$REF_GFF" out="$POSITIVES_GFF" op=intersect

# --- 3.3: Convert to TSV and label ---
echo "[3/4] Converting GFF sets to labeled TSV format..."
NEGATIVES_TSV="$TMPDIR/negatives.tsv"
POSITIVES_TSV="$TMPDIR/positives.tsv"

# Convert negatives with label -1
gff2tsv.sh in="$NEGATIVES_GFF" out="$NEGATIVES_TSV" label=-1
# Convert positives with label 1
gff2tsv.sh in="$POSITIVES_GFF" out="$POSITIVES_TSV" label=1

# --- 3.4: Combine into final training set ---
echo "[4/4] Assembling final training set: $OUT_TSV"
# Create the final file and add a header (optional but good practice)
echo "contig\tlog_len\t...\tmatch" > "$OUT_TSV"
cat "$POSITIVES_TSV" >> "$OUT_TSV"
cat "$NEGATIVES_TSV" >> "$OUT_TSV"

# --- Cleanup ---
rm -r "$TMPDIR"
echo "Pipeline complete."
