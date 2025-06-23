#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_6dims_400k_ncbi_100_6f.log
#SBATCH --error=cycles_trainer_6dims_400k_ncbi_100_6f.err
#SBATCH --job-name=cycles_trainer_6dims_400k_ncbi_100_6f

echo "Starting cycles training at $(date)"

./bbmap/train.sh in=ncbi_100_6f.tsv out=cycles_trained_6dims_400k_ncbi_100_6f.bbnet maxdims=354,300,240,180,90,1 mindims=354,150,120,90,10,1 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
