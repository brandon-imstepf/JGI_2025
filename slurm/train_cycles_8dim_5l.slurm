#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_100_8f_5l.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_100_8f_5l.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_100_8f_5l

echo "Starting cycles training at $(date)"

../bbmap/train.sh in=../ncbi_100_8f.tsv out=cycles_trained_8dims_400k_ncbi_100_8f_5l.bbnet maxdims=356,300,240,180,1 mindims=356,150,120,90,1 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
