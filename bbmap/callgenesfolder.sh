#!/bin/bash

# Require in= and out= as first two arguments
if [[ "$1" != in=* || "$2" != out=* ]]; then
    echo "Usage: $0 in=<genefolder> out=<outvector> [additional CallGenes args, separated by spaces]"
    exit 1
fi

GENEFOLDER="${1#in=}"
OUTVECTOR="${2#out=}"
shift 2
ARGS="$@"

if [[ -z "$GENEFOLDER" || -z "$OUTVECTOR" ]]; then
    echo "Usage: $0 in=<genefolder> out=<outvector> [additional CallGenes args]"
    exit 1
fi

TMPDIR=$(mktemp -d)
echo "Temporary directory: $TMPDIR"

# Step 1: Run CallGenes for each pair
for fna in "$GENEFOLDER"/*.fna.gz; do
    base=$(basename "$fna" .fna.gz)
    gff="$GENEFOLDER/$base.gff.gz"
    if [[ -f "$gff" ]]; then
        echo "Processing $base with ml, cds, $ARGS"
        # Output to a temp file for each pair
        java -cp /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/bbmap/current/ prok.CallGenes in="$fna" truegenes="$gff" outvector="$TMPDIR/$base.outvector" ml cds $ARGS -Xmx1g 
    else
        echo "WARNING: No matching .gff.gz for $fna, skipping."
    fi
done

# Step 2: Concatenate all outvector files, keeping only the first header
first=1
for file in "$TMPDIR"/*.outvector; do
    if [[ $first -eq 1 ]]; then
        cat "$file" > "$OUTVECTOR"
        first=0
    else
        # Skip header (assumes header is the first line)
        tail -n +2 "$file" >> "$OUTVECTOR"
    fi
done

# Step 3: Clean up
rm -r "$TMPDIR"
echo "All done. Output in $OUTVECTOR"
