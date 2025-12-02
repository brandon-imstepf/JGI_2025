#!/bin/bash

usage(){
echo "
Written by Brandon Imstepf
Last modified July 21st, 2025

Description:  Trains or evaluates a convolutional neural network for binary classification.

Usage:  cnntrain.sh in=<data> out=<network.bbnet> filters=32,64 dense=128 epochs=10

I/O parameters:
in=<file>       Tab-delimited data vectors. The first line should look like
data=<file>     '#dims 356 1' with the number of inputs and outputs; the
                first X columns are inputs, and the last Y the desired result.
                Subsequent lines are tab-delimited floating point numbers.
validate=<file> Optional validation dataset used exclusively for evaluation.
                If not provided, 10% of training data is used for validation.
net=<file>      Optional input network for evaluation mode or continued training.
out=<file>      Final output network after training (*.bbnet format).
maxlines=<int>  Maximum training samples to load (default: all).
maxlinesv=<int> Maximum validation samples to load (default: all).

Architecture parameters:
filters=X,Y,Z   Number of filters in each convolutional layer (default: 32,64).
                Example: filters=32,64,128 creates 3 conv layers.
filtersize=X,Y  Size of convolutional filters (default: 5,3).
                If fewer sizes than layers, last size is repeated.
poolsize=X,Y    Size of max pooling windows (default: 2,2).
                Pooling occurs after each conv layer.
dense=X,Y       Number of neurons in dense layers (default: 128).
                Example: dense=256,128 creates 2 dense layers.

Training parameters:
epochs=<int>    Number of training epochs (default: 10).
batch=<int>     Batch size for mini-batch training (default: 32).
lr=<float>      Learning rate (default: 0.001).
dropout=<float> Dropout rate for regularization (default: 0.5).
shuffle         Shuffle training data before splitting.
split=<float>   Fraction to use for validation if validate not specified (default: 0.1).

Evaluation mode:
evaluate        Run in evaluation mode instead of training.
                Requires net=<file> with trained network.
                Outputs comprehensive metrics including ROC analysis.

Output control:
verbose         Show per-batch training progress.
out=<file>      Write console output to file.
                Use 'stdout' or 'stderr' for standard streams.

Example usage:
# Train a new network
cnntrain.sh data=train.tsv out=model.bbnet filters=32,64,128 dense=256,128 epochs=20

# Evaluate on test data
cnntrain.sh data=test.tsv net=model.bbnet evaluate

# Train with custom architecture
cnntrain.sh data=train.tsv validate=val.tsv out=model.bbnet \\
    filters=16,32,64,128 filtersize=7,5,3,3 poolsize=2,2,2,2 \\
    dense=512,256,128 epochs=50 batch=64 lr=0.0001

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.
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
