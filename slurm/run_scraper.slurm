#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 24:00:00            # Time
#SBATCH -c 1                   # 1 CPU core (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory
#SBATCH --output=ncbi_scraper.log
#SBATCH --error=ncbi_scraper.err
#SBATCH --job-name=ncbi_scraper

# Load your Python environment
module load python
source /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/jgi_venv/bin/activate

python3 /clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/NCBI_Scraper.py
