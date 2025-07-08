#!/bin/bash
# A simple wrapper for the GffToTsv.java
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
z="-Xmx1g"

java $z -cp $CP prok.GffToTsv "$@"
