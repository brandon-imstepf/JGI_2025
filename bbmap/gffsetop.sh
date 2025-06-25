#!/bin/bash

usage(){
echo "
Written by Brandon Imstepf
Last modified June 25, 2025

Description: 
Performs set operations on the features of two GFF files.
It identifies features as identical if their contig, start, stop,
and strand are the same. All other fields are ignored for comparison
but preserved from the first file (in_a) for the output.

Usage:   gffsetop.sh in_a=<file1.gff> in_b=<file2.gff> out=<outfile.gff> op=<operation>

Required Parameters:
in_a=<file>     The first GFF file (e.g., the set to subtract from).
in_b=<file>     The second GFF file.
out=<file>      The resulting output GFF file.
op=<string>     The set operation to perform. Options are:
                - intersect:  Features present in BOTH file A and file B.
                - union:      All unique features from file A and file B combined.
                - subtract:   Features present in file A but NOT in file B.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

# Assumes 'current' directory is in the same dir as the script
CP="$DIR""current/"

# This block calculates memory settings automatically.
z="-Xmx4g"
z2="-Xms4g"
set=0
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
    usage
    exit
fi
calcXmx () {
    source "$DIR""/calcmem.sh"
    setEnvironment
    parseXmx "$@"
}
calcXmx "$@"

function gffsetop() {
    # We're calling a new Java class: prok.GffSetOperator
    local CMD="java $EA $EOOM $z $z2 -cp $CP prok.GffSetOperator $@"
    echo "Executing command:" >&2
    echo "$CMD" >&2
    eval $CMD
}

gffsetop "$@"
