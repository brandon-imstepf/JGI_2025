#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 24:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_5dims_400k_d075_ncbi_100_6f.log
#SBATCH --error=cycles_trainer_5dims_400k_d075_ncbi_100_6f.err
#SBATCH --job-name=cycles_trainer_5dims_400k_d075_ncbi_100_6f

echo "Starting cycles training at $(date)"

./bbmap/train.sh in=ncbi_100_6f.tsv out=cycles_trained_5dims_400k_d075_ncbi_100_6f.bbnet maxdims=353,270,180,90,1 mindims=353,100,90,10,1 cycles=4 nets=8 density=0.75 -Xmx400g

echo "Finished cycles training at $(date)"
