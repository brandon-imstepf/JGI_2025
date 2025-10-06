#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cnntrain_filters-32-64-128_dense-256-128_epochs-10_7-23.log
#SBATCH --error=cnntrain_filters-32-64-128_dense-256-128_epochs-10_7-23.err
#SBATCH --job-name=cnntrain_filters-32-64-128_dense-256-128_epochs-10

echo "Starting CNN training at $(date)"

../bbmap/cnntrain.sh in=../ncbi_10_balanced.tsv in=../ncbi_10_balanced.tsv out=cnn_filters-32-64-128_dense-256-128_epochs-10.bbnet filters=32,64,128 dense=256,128 epochs=10

echo "Finished CNN training at $(date)"
