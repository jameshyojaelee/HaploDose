#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --job-name=masld-test
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Runs submit_cellranger.sh every 30 min until every sample in SAMPLES_TSV is submitted (it submits a patient's
# samples only when all of that patient's FASTQs have .ok). Submits nothing else.
# Usage: sbatch --output=OUT/logs/watch_%j.out watch_cellranger.sbatch OUT SAMPLES_TSV
set -euo pipefail
OUT=$1 TABLE=$2
total=$(($(wc -l < "$TABLE") - 1))
while true; do
  date -Is
  bash benchmarks/spectrum/submit_cellranger.sh "$OUT" "$TABLE"
  (( $(wc -l < "$OUT/jobs.tsv") - 1 >= total )) && break
  sleep 1800
done
