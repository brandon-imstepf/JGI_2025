#!/bin/bash

usage(){
echo "
Written by  Brandon Imstepf
Last modified July 7th, 2025

Description:  Trains a convolutional neural network.

Usage:  cnntrain.sh in=<data>  out=<trained network>

I/O parameters:
in=<file>       Tab-delimited data vectors.  The first line should look like
                '#dims	5	1' with the number of inputs and outputs; the
                first X columns are inputs, and the last Y the desired result.
                Subsequent lines are tab-delimited floating point numbers.

validate=<file> Optional validation dataset used exclusively for evaluation.
net=<file>      Optional input network to train.
out=<file>      Final output network after the last epoch.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

z="-Xmx8g"
z2="-Xms8g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 8000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

train() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP ml.CNNTrainer $@"
	echo $CMD >&2
	eval $CMD
}

train "$@"
