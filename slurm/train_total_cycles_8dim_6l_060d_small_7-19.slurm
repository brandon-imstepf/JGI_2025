#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 36:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_100_8f_6l_060d_small_7-19.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_100_8f_6l_060d_small_7-19.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_100_8f_6l_060d_small_7-19

echo "Starting cycles training at $(date)"

../bbmap/train.sh in=../ncbi_total_training_final.tsv out=cycles_trained_8dims_400k_ncbi_100_8f_6l_060d_small_7-19.bbnet maxdims=356,150,200,90,12,1 mindims=356,50,100,40,5,1 density=0.60 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
