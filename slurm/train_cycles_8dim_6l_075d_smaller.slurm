#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_100_8f_6l_075d_smaller.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_100_8f_6l_075d_smaller.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_100_8f_6l_075d_smaller

echo "Starting cycles training at $(date)"

../bbmap/train.sh in=../ncbi_100_8f.tsv out=cycles_trained_8dims_400k_ncbi_100_8f_6l_075d_smaller.bbnet maxdims=356,100,250,75,10,1 mindims=356,50,100,40,5,1 density=0.75 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
