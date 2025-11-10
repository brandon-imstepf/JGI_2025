#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH -c 8                   # 8 CPU cores (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory
#SBATCH --output=ncbi_builder_11-7.log
#SBATCH --error=ncbi_builder_11-7.err
#SBATCH --job-name=ncbi_builder


echo "Started building at $(date)"

btsf_new.sh in=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/NCBI_Datasets_Training out=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/ncbi_master_11-7.tsv

echo "Finished building at $(date)"