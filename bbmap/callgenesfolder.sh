#!/bin/bash

usage(){
echo "
Written by Brandon Imstepf, inspired by Brian Bushnell's scripts.
Last modified June 24, 2025

Description: 
Processes a folder of FASTA files to generate gene-call data.
It can operate in two modes depending on whether the 'ml' flag is passed.
- GFF Mode (default): Creates an annotated GFF file for each input FASTA,
  then concatenates them into a single GFF output file.
- ML Mode ('ml' flag): Creates a feature vector file for each input,
  then concatenates them into a single TSV output file.

Usage:   callgenesfolder.sh in=<folder> out=<outfile> [options]

Required Parameters:
in=<folder>     Input folder containing .fna.gz files.
                For each sample.fna.gz, a sample.gff.gz is expected if using 'truegenes'.
out=<outfile>   Final, combined output file (.gff or .tsv).

Common Pass-through Options (for prok.CallGenes):
ml              Activates Machine Learning mode, outputting a feature vector TSV.
cds             Call protein-coding genes.
truegenes=      Specify a reference GFF for labeling training data.
                Note: this script automatically finds the matching GFF for each FASTA.

Example (GFF Mode):
callgenesfolder.sh in=./genomes out=all_annotations.gff cds

Example (ML Mode):
callgenesfolder.sh in=./genomes out=all_vectors.tsv ml cds
"
}


# This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

CP="$DIR""current/" # Assumes 'current' directory is in the same dir as the script

# This block calculates memory settings automatically.
z="-Xmx6g"
z2="-Xms6g"
set=0
calcXmx () {
    source "$DIR""/calcmem.sh"
    setEnvironment
    parseXmx "$@"
}
calcXmx "$@"


# --- Argument Parsing ---
INFOLDER=""
OUTFILE=""
PASSTHRU_ARGS=() # Use an array to safely handle arguments with spaces

for arg in "$@"; do
    case $arg in
        in=*)
            INFOLDER="${arg#in=}"
            ;;
        out=*)
            OUTFILE="${arg#out=}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            # Anything else is a pass-through argument for the Java program
            PASSTHRU_ARGS+=("$arg")
            ;;
    esac
done

# Check if required arguments were provided
if [[ -z "$INFOLDER" || -z "$OUTFILE" ]]; then
    echo "ERROR: Both in=<folder> and out=<outfile> are required." >&2
    usage
    exit 1
fi


# --- Main Logic ---

# Check if 'ml' is in the pass-through arguments to determine the mode
ML_MODE=0
for arg in "${PASSTHRU_ARGS[@]}"; do
    if [[ "$arg" == "ml" ]]; then
        ML_MODE=1
        break
    fi
done

TMPDIR=$(mktemp -d)
echo "Temporary directory: $TMPDIR"

# Step 1: Run CallGenes for each FASTA file
for fna in "$INFOLDER"/*.fna.gz; do
    base=$(basename "$fna" .fna.gz)
    # Define a temporary output file for this specific run
    TEMP_OUT="$TMPDIR/$base.tmp"

    # Assemble the command for this specific file
    CMD_ARGS=("in=$fna")
    
    # Check for a matching GFF file if truegenes might be used
    gff="$INFOLDER/$base.gff.gz"
    if [[ -f "$gff" ]]; then
        CMD_ARGS+=("truegenes=$gff")
    fi

    if [[ $ML_MODE -eq 1 ]]; then
        echo "Processing (ML Mode): $base"
        CMD_ARGS+=("outvector=$TEMP_OUT")
    else
        echo "Processing (GFF Mode): $base"
        CMD_ARGS+=("outgff=$TEMP_OUT")
    fi
    
    # Run the Java application with the combined arguments
    java $EA $EOOM $z $z2 -cp $CP prok.CallGenes "${CMD_ARGS[@]}" "${PASSTHRU_ARGS[@]}"
done


# Step 2: Concatenate all temporary output files into the final destination
echo "Combining results into $OUTFILE"
# Clear the output file to start fresh
> "$OUTFILE" 
# Use find to handle cases with no temp files gracefully
find "$TMPDIR" -name '*.tmp' -print0 | xargs -0 cat >> "$OUTFILE"


# Step 3: Clean up
rm -r "$TMPDIR"
echo "All done."