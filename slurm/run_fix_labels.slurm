#!/bin/bash
#SBATCH -A grp-org-gt-ga
#SBATCH -q jgi_normal
#SBATCH -t 24:00:00
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --job-name=fix_labels
#SBATCH --output=fix_labels.%j.log
#SBATCH --error=fix_labels.%j.err

set -euo pipefail

ROOT="/clusterfs/jgi/scratch/gentech/genome_analysis/brandonimstepf/training_datasets"
AWK_PROG='BEGIN{OFS="\t"}
{
    last=$NF
    gsub(/\x00/, "", last)
    if (last ~ /^1+$/) { $NF = 1 }
    else if (last ~ /^0+$/) { $NF = 0 }
    else { $NF = substr(last,1,1) }
    print
}'

if [ "$#" -gt 0 ]; then
  dirs=("$@")
else
  dirs=("$ROOT"/*)
fi

for dir in "${dirs[@]}"; do
  [ -d "$dir" ] || continue
  echo "[start] $dir"
  date
  shopt -s nullglob
  files=("$dir"/*.tsv.gz)
  shopt -u nullglob
  if [ ${#files[@]} -eq 0 ]; then
    echo "  no TSVs in $dir"
    echo "[done ] $dir"
    date
    continue
  fi
  for file in "${files[@]}"; do
    echo "  processing $file"
    tmp="${file%.gz}.tmp.gz"
    zcat "$file" | awk "$AWK_PROG" | pigz -c > "$tmp"
    mv "$tmp" "$file"
  done
  echo "[done ] $dir"
  date
  echo
  if [ $SECONDS -ge $((24*3600-600)) ]; then
    echo "Reached time budget (~24h); exiting."
    exit 0
  fi
  sleep 5
 done
