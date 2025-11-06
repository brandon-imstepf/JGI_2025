#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH -c 8                   # 8 CPU cores (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory
#SBATCH --output=ncbi_builder_10-23_orig.log
#SBATCH --error=ncbi_builder_10-23_orig.err
#SBATCH --job-name=ncbi_builder_orig


echo "Started building at $(date)"

btsf_new.sh in=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/NCBI_Datasets_Training out=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/ncbi_master_11-5_orig.tsv

echo "Finished building at $(date)"