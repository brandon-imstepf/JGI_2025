#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_7-21.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_7-21.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_7-21

echo "Starting cycles training at $(date)"

../bbmap/train.sh in=../ncbi_total_training_final.tsv out=cycles_trained_8dims_400k_ncbi_total_8f_6l_small_7-21.bbnet maxdims=356,125,275,75,25,1 mindims=356,50,150,25,5,1 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
