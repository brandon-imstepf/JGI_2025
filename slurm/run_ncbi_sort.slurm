#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH -c 8                   # 8 CPU cores (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory
#SBATCH --output=ncbi_sorter_10-29.log
#SBATCH --error=ncbi_sorter_10-29.err
#SBATCH --job-name=ncbi_sorter

echo "Starting sorting at $(date)"

/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/bbmap/buildtrainingset_ncbi80_10_10_v2.sh in=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/NCBI_Datasets_10-23_SORTED out=/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/10-29_master_data.tsv seed=42

echo "Finished sorting at $(date)"