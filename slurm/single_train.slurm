#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 24:00:00            # Time
#SBATCH -c 8                   # 1 CPU core (adjust if you want more)
#SBATCH --mem=16G              # 16 GB memory
#SBATCH --output=single_trainer.log
#SBATCH --error=single_trainer.err
#SBATCH --job-name=single_trainer


./bbmap/train.sh in=ncbi_vector_long.tsv dims=353,300,200,100,1 outb=353-300-200-100.bbnet

