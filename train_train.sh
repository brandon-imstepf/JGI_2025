#!/bin/bash

# Input file
input="ncbi_vector_long.tsv"

# Arrays of dims and out parameters

dims_list=("353,180,90,45,1" "353,45,90,180,1" "353,90,45,180,1" "353,300,250,200,150,100,50,1" "353,100,50,200,1")
out_list=("out_353-180-90-45-1.bbnet" "out_353-45-90-180-1.bbnet" "out_353-90-45-180-1.bbnet" "out_353-300-250-200-150-100-50-1.bbnet" "out_353-100-50-200-1.bbnet")

for i in "${!dims_list[@]}"; do
  dims="${dims_list[$i]}"
  out="${out_list[$i]}"
  jobname="train_${dims//,/}" # e.g. train_35390451801
  outlog="${out}.out.log"
  errlog="${out}.err.log"
  export DIMS="${dims}"
  sbatch --job-name="$jobname" \
         --output="$outlog" \
         --error="$errlog" \
         --export=ALL,INPUT="$input",DIMS,OUT="$out" \
         train_train.slurm
  echo "Submitted job: $jobname with dims=$dims, outb=$out"
done
