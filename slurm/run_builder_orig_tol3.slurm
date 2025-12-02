#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --output=ncbi_builder_11-7_orig_tol3.log
#SBATCH --error=ncbi_builder_11-7_orig_tol3.err
#SBATCH --job-name=ncbi_builder_orig_tol3

echo "Started building (orig tol=3) at $(date)"

btsf_new.sh \
    in=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/NCBI_Datasets_Training \
    out=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/ncbi_master_11-7_orig_tol3.tsv \
    label_start_tolerance=3

echo "Finished building (orig tol=3) at $(date)"
