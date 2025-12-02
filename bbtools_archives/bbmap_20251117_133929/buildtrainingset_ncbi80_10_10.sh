#!/bin/bash
# buildtrainingset_ncbi80_10_10.sh
#
# Created to implement the revised NN training strategy (80/10/10 split,
# strict no-overlap negatives against the global NCBI gene set).
#
# For each genome (pair of .fna.gz and .gff.gz) in the input folder the
# script:
#  1) Calls `callgenes` in nofilter mode to produce ALL candidate ORFs.
#  2) Calls `callgenes` in normal (filtered) mode to produce high-confidence
#     predictions (used to identify true positives where they intersect
#     with NCBI annotations).
#  3) Concatenates NCBI GFFs into a global `ncbi_all.gff` and builds a
#     global mask of NCBI gene intervals.
#  4) Concatenates ALL candidate ORFs into `all_orfs_all.gff` and removes
#     any ORF that overlaps ANY NCBI gene (global) to obtain negatives.
#  5) Produces a deterministic, seeded contig->split mapping (80/10/10 by
#     cumulative basepairs) and splits positives/negatives by contig.
#  6) Converts split GFF files to TSV (labeled) and balances negatives
#     (max 2x positives per split by default).
#
# Usage: buildtrainingset_ncbi80_10_10.sh in=<folder> out=<master.tsv> [options]
# Options:
#   seed=<int>            Deterministic random seed (default: 42)
#   split_by=contig|genome  Unit to split by; default: contig (recommended)
#   max_neg_ratio=<float> Max negatives per positives (default: 2.0)
#   keep_tmp=<bool>       Keep temporary working directory (default false)
#
# Notes:
# - This script expects the following helper scripts to be available in the
#   same repo directory: callgenes.sh, gffsetop.sh, gff2tsv.sh
# - If `bedtools` is available it will be used for interval subtraction; if
#   not available this script falls back to `gffsetop.sh` which must support
#   op=subtract between GFF files.

set -euo pipefail

usage(){
  cat <<'USAGE'
Usage: buildtrainingset_ncbi80_10_10.sh in=<folder> out=<master.tsv> [options]

Required:
  in=<folder>       Input folder containing matching .fna.gz and .gff.gz files
  out=<master.tsv>  Output combined TSV training file

Optional:
  seed=42           Deterministic seed (default 42)
  split_by=contig   Unit to split by: contig or genome (default contig)
  max_neg_ratio=2.0 Max allowed negatives per positive (default 2.0)
  keep_tmp=false    Keep temporary directory for debugging (default false)
  threads=8         Number of threads passed to callgenes (default 8)
\nExamples:
  buildtrainingset_ncbi80_10_10.sh in=./data out=master.tsv seed=123
USAGE
}

# --- Boilerplate to find repo scripts (same pattern as other scripts) ---
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

# --- Defaults ---
INFOLDER=""
OUT_TSV=""
SEED=42
SPLIT_BY="contig"
MAX_NEG_RATIO=2.0
KEEP_TMP=false
THREADS=8

# --- Parse args ---
PASSTHRU_ARGS=()
for arg in "$@"; do
  case $arg in
    in=*) INFOLDER="${arg#in=}";;
    out=*) OUT_TSV="${arg#out=}";;
    seed=*) SEED="${arg#seed=}";;
    split_by=*) SPLIT_BY="${arg#split_by=}";;
    max_neg_ratio=*) MAX_NEG_RATIO="${arg#max_neg_ratio=}";;
    keep_tmp=*) KEEP_TMP="${arg#keep_tmp=}";;
    threads=*) THREADS="${arg#threads=}";;
    -h|--help) usage; exit 0;;
    *) PASSTHRU_ARGS+=("$arg");;
  esac
done

if [[ -z "$INFOLDER" || -z "$OUT_TSV" ]]; then
  echo "ERROR: Missing required parameters in= and out=" >&2
  usage
  exit 1
fi

# --- Temporary workspace ---
TMPDIR=$(mktemp -d -t buildtrainingset.XXXX)
echo "Using temporary directory: $TMPDIR"
trap 'if [ "$KEEP_TMP" = "false" ]; then rm -rf "$TMPDIR"; fi' EXIT

ALL_ORFS_GLOBAL="$TMPDIR/all_orfs_all.gff"
NCBI_ALL_GFF="$TMPDIR/ncbi_all.gff"
NCBI_BED="$TMPDIR/ncbi_intervals.bed"
CONTIG_LENGTHS="$TMPDIR/contig_lengths.tsv"
CONTIG_SPLIT="$TMPDIR/contig_to_split.csv"

mkdir -p "$TMPDIR/per_genome"

echo "Gathering inputs from $INFOLDER"

# Concatenate all NCBI GFFs (assume files named <base>.gff.gz)
rm -f "$NCBI_ALL_GFF" "$ALL_ORFS_GLOBAL"

for ref_fasta in "$INFOLDER"/*.fna.gz; do
  base=$(basename "$ref_fasta" .fna.gz)
  ref_gff="$INFOLDER/$base.gff.gz"
  if [[ ! -f "$ref_gff" ]]; then
    echo "WARNING: No matching reference GFF for $ref_fasta. Skipping."
    continue
  fi

  echo "--- Processing genome: $base ---"

  # Per-genome temporary files
  ALL_GFF="$TMPDIR/per_genome/${base}.all.gff"
  PRED_GFF="$TMPDIR/per_genome/${base}.predicted.gff"

  # 1) Generate ALL candidates (nofilter=true) and filtered predictions
  echo "  Calling genes (all candidates)..."
  "$DIR/callgenes.sh" in="$ref_fasta" outgff="$ALL_GFF" cds nofilter=true threads=$THREADS "${PASSTHRU_ARGS[@]}"

  echo "  Calling genes (filtered predictions)..."
  "$DIR/callgenes.sh" in="$ref_fasta" outgff="$PRED_GFF" cds threads=$THREADS "${PASSTHRU_ARGS[@]}"

  # Append per-genome ALL candidates to global file
  if [[ -s "$ALL_GFF" ]]; then
    cat "$ALL_GFF" >> "$ALL_ORFS_GLOBAL"
  fi

  # Append NCBI reference for global NCBI file (decompress if needed)
  if [[ "$ref_gff" == *.gz ]]; then
    zcat "$ref_gff" >> "$NCBI_ALL_GFF"
  else
    cat "$ref_gff" >> "$NCBI_ALL_GFF"
  fi

  # Create contig lengths for this genome from fasta
  # Produce lines: contig<TAB>length<TAB>genome_base
  # Normalize the contig id by trimming at first whitespace to match GFF contig IDs
  zcat "$ref_fasta" | awk -v base="$base" 'BEGIN{RS=">";FS="\n"} NR>1{sum=0; for(i=2;i<=NF;i++) sum+=length($i); split($1,a,/\s+/); id=a[1]; print id "\t" sum "\t" base}' >> "$CONTIG_LENGTHS"

done

if [[ ! -s "$NCBI_ALL_GFF" ]]; then
  echo "ERROR: No NCBI GFF files were found in $INFOLDER" >&2
  exit 1
fi

echo "Building NCBI interval mask (BED)"
# Convert NCBI GFF -> BED (0-based, half-open). Keep only gene/CDS lines
awk '$0!~/^#/ && NF>8 {print $1"\t"($4-1)"\t"$5"\t"$3}' "$NCBI_ALL_GFF" > "$NCBI_BED"

echo "Filtering all ORFs against global NCBI gene set to make negatives"
# Use the repository tool `gffsetop.sh` to subtract NCBI genes from the full ORF set.
# This avoids any dependency on external bedtools.
touch "$ALL_ORFS_GLOBAL"
"$DIR/gffsetop.sh" in_a="$ALL_ORFS_GLOBAL" in_b="$NCBI_ALL_GFF" out="$TMPDIR/all_orfs_no_ncbi_overlap.gff" op=subtract

echo "Creating deterministic contig->split mapping (seed=$SEED, split_by=$SPLIT_BY)"
# Use the contig lengths file to create a weighted split by cumulative bp
python3 - <<PY
import csv,random
random.seed(int($SEED))
lines=[]
with open('$CONTIG_LENGTHS') as f:
    for r in csv.reader(f, delimiter='\t'):
        if len(r)<2: continue
        contig=r[0]; length=int(r[1]); genome=r[2] if len(r)>2 else ''
        lines.append((contig,length,genome))
if not lines:
    raise SystemExit('No contigs found in $CONTIG_LENGTHS')
total=sum(l for _,l,_ in lines)
random.shuffle(lines)
cuts=(0.8,0.9,1.0)
cum=0
out=[('contig','split','genome')]
for contig,length,genome in lines:
    frac=(cum+length)/total
    if frac<=cuts[0]: s='train'
    elif frac<=cuts[1]: s='val'
    else: s='test'
    out.append((contig,s,genome))
    cum+=length
with open('$CONTIG_SPLIT','w') as fw:
    for row in out:
        fw.write(','.join(row)+'\n')
print('WROTE $CONTIG_SPLIT')
PY

echo "Splitting NCBI positives and negatives by contig mapping"
# Prepare per-split files
NCBI_TRAIN="$TMPDIR/ncbi_all_train.gff"
NCBI_VAL="$TMPDIR/ncbi_all_val.gff"
NCBI_TEST="$TMPDIR/ncbi_all_test.gff"
NEG_TRAIN="$TMPDIR/all_orfs_no_ncbi_overlap_train.gff"
NEG_VAL="$TMPDIR/all_orfs_no_ncbi_overlap_val.gff"
NEG_TEST="$TMPDIR/all_orfs_no_ncbi_overlap_test.gff"
: > "$NCBI_TRAIN"; : > "$NCBI_VAL"; : > "$NCBI_TEST"
: > "$NEG_TRAIN"; : > "$NEG_VAL"; : > "$NEG_TEST"

python3 - <<PY
import csv
split={}
with open('$CONTIG_SPLIT') as f:
    r=csv.DictReader(f)
    for row in r:
        split[row['contig']]=row['split']
def route_gff(infile, outprefix):
    outs={'train':open(outprefix+'_train.gff','a'),'val':open(outprefix+'_val.gff','a'),'test':open(outprefix+'_test.gff','a')}
    for line in open(infile):
        if line.startswith('#') or not line.strip(): continue
        contig=line.split()[0]
        s=split.get(contig)
        if s:
            outs[s].write(line)
for infile,outprefix in [('$NCBI_ALL_GFF','$NCBI_TRAIN'.replace('_train.gff','') ), ('$TMPDIR/all_orfs_no_ncbi_overlap.gff','$NEG_TRAIN'.replace('_train.gff',''))]:
    route_gff(infile,outprefix)
print('Split GFFs written')
PY

echo "Converting GFFs to TSV (labeling: positives=1, negatives=0)"
# Convert positives
"$DIR/gff2tsv.sh" in="$NCBI_TRAIN" out="$TMPDIR/ncbi_train.tsv" label=1 || true
"$DIR/gff2tsv.sh" in="$NCBI_VAL" out="$TMPDIR/ncbi_val.tsv" label=1 || true
"$DIR/gff2tsv.sh" in="$NCBI_TEST" out="$TMPDIR/ncbi_test.tsv" label=1 || true
# Convert negatives
"$DIR/gff2tsv.sh" in="$NEG_TRAIN" out="$TMPDIR/neg_train.tsv" label=0 || true
"$DIR/gff2tsv.sh" in="$NEG_VAL" out="$TMPDIR/neg_val.tsv" label=0 || true
"$DIR/gff2tsv.sh" in="$NEG_TEST" out="$TMPDIR/neg_test.tsv" label=0 || true

echo "Balancing negatives per split (max ratio = $MAX_NEG_RATIO)"
for split in train val test; do
  POS="$TMPDIR/ncbi_${split}.tsv"
  NEG="$TMPDIR/neg_${split}.tsv"
  if [[ ! -f "$POS" ]]; then
    POS_COUNT=0
  else
    POS_COUNT=$(wc -l < "$POS" || echo 0)
  fi
  if [[ ! -f "$NEG" ]]; then
    NEG_COUNT=0
  else
    NEG_COUNT=$(wc -l < "$NEG" || echo 0)
  fi
  MAX_NEG=$(printf "%.0f" $(echo "$POS_COUNT * $MAX_NEG_RATIO" | bc -l))
  if [[ $NEG_COUNT -gt $MAX_NEG && $MAX_NEG -gt 0 ]]; then
    echo "  Subsampling negatives for $split: $NEG_COUNT -> $MAX_NEG"
    shuf -n "$MAX_NEG" "$NEG" > "$TMPDIR/neg_${split}.balanced.tsv"
    mv "$TMPDIR/neg_${split}.balanced.tsv" "$NEG"
  fi
done

echo "Assembling master TSV: $OUT_TSV"
echo "# Combined training set produced by buildtrainingset_ncbi80_10_10.sh" > "$OUT_TSV"
# Append train/val/test in that order (optional: could shuffle later)
for split in train val test; do
  for f in "$TMPDIR/ncbi_${split}.tsv" "$TMPDIR/neg_${split}.tsv"; do
    if [[ -f "$f" ]]; then
      cat "$f" >> "$OUT_TSV"
    fi
  done
done

echo "Done. Master TSV: $OUT_TSV"
if [ "$KEEP_TMP" = "true" ]; then
  echo "Temporary files retained at $TMPDIR"
else
  echo "Temporary files removed"
fi
