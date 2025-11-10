#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_orig

echo "Starting cycles training at $(date)"

../bbmap/train.sh in=../ncbi_master_11-7_splits/ncbi_master_11-7_orig_train.tsv out=cycles_trained_8dims_400k_ncbi_total_8f_6l_small_11-7_orig.bbnet maxdims=358,125,275,75,25,1 mindims=358,50,150,25,5,1 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
