#!/bin/bash

#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 24:00:00            # Time
#SBATCH -c 8                   # 1 CPU core (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory


echo "Running on input: $INPUT"
echo "Dims: $DIMS"
echo "Output: $OUT"
./bbmap/train.sh in="$INPUT" dims="$DIMS" outb="$OUT" balance=f t=8

