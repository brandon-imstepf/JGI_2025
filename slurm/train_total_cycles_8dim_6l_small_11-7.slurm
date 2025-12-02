#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 48:00:00            # Time
#SBATCH --exclusive
#SBATCH --output=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_new.log
#SBATCH --error=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_new.err
#SBATCH --job-name=cycles_trainer_8dims_400k_ncbi_total_8f_6l_small_11-7_new

echo "Starting cycles training at $(date)"

INPUT_TSV=../ncbi_master_11-7_splits/ncbi_master_11-7_train.tsv
SHUFFLED_TSV=../ncbi_master_11-7_splits/ncbi_master_11-7_train.shuffled.tsv

echo "Shuffling training data into ${SHUFFLED_TSV}"
python3 ../python/shuffle_tsv.py --input "${INPUT_TSV}" --output "${SHUFFLED_TSV}" --chunk 500000

../bbmap/train.sh in="${SHUFFLED_TSV}" out=cycles_trained_8dims_400k_ncbi_total_8f_6l_small_11-7_new.bbnet maxdims=358,125,275,75,25,1 mindims=358,50,150,25,5,1 cycles=4 nets=8 -Xmx400g

echo "Finished cycles training at $(date)"
