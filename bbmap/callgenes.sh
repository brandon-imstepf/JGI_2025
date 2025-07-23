#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 23, 2025

Description:  Finds orfs and calls genes in unspliced prokaryotes.
This includes bacteria, archaea, viruses, and mitochondria.
Can also predict 16S, 18S, 23S, 5S, and tRNAs.
Includes an optional neural network filtering step for improved accuracy.

Usage:  callgenes.sh in=contigs.fa out=calls.gff
NN Filtering:  callgenes.sh in=contigs.fa out=filtered.gff net=model.bbnet

File parameters:
in=<file>       A fasta file; the only required parameter.
out=<file>      Output gff file.
outa=<file>     Amino acid output.
out16s=<file>   16S output.
model=<file>    A pgm file or comma-delimited list.
               If unspecified a default model will be used.
stats=stderr    Stats output (may be stderr, stdin, a file, or null).
hist=null       Gene length histogram.
compareto=      Optional reference gff file to compare with the gene calls.
               'auto' will name it based on the input file name.

Neural Network Filtering Parameters:
net=<file>      Use a neural network to filter CDS candidates. When this flag
               is used, the program runs in inference mode.
cutoff=0.5     Score cutoff for keeping a gene candidate.
lowpass=f      Set to true to keep genes with scores *below* the cutoff.

Training Data Generation Parameters:
nofilter=f      Generate all possible gene candidates, bypassing internal and
               NN filters. Used for creating training sets.
truegenes=<file> Provide a GFF of known genes to label candidates as true (1)
               or false (0) in the output vector.
seq=f           Output the raw DNA sequence of candidates instead of one-hot
               encoded windows.

Formatting parameters:
json=false      Print stats in JSON.
binlen=21       Histogram bin length.
bins=1000      Maximum histogram bins.
pz=f           (printzero) Print histogram lines with zero count.

Other parameters:
minlen=60       Don't call genes shorter than this.
trd=f           (trimreaddescription) Set to true to trim read headers after
               the first whitespace.
merge=f        For paired reads, merge before calling.
detranslate=f  Output canonical nucleotide sequences instead of amino acids.
recode=f       Re-encode nucleotide sequences over called genes, leaving
               non-coding regions unchanged.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx6g"
z2="-Xms6g"
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

function callgenes() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP prok.CallGenes $@"
	#Too long to echo sometimes since wildcards can be expanded
	#echo $CMD >&2
	eval $CMD
}

callgenes "$@"
