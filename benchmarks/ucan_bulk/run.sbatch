#!/usr/bin/env bash
#SBATCH --job-name=ucan_bulk
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_bulk_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_bulk_%j.err
# U-CAN bulk RNA benchmark: Haplodose bulk, naive log ratio, infercnvpy; scored against WGS FACETS.
# Usage: sbatch run.sbatch INPUT_DIR ROOT
set -euo pipefail
IN=$1 ROOT=$2
mkdir -p "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
B=benchmarks/ucan_bulk
test -f "$ROOT/haplodose/run.json" || /usr/bin/time -v $HD -s -m haplodose.cli bulk --counts "$IN/counts.tsv" \
  --samples "$IN/samples.tsv" --positions "$IN/positions.tsv" \
  --exclude-genes /gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt --out "$ROOT/haplodose" 2> "$ROOT/haplodose_time.txt"
test -f "$ROOT/naive_arms.tsv" || $HD -s $B/naive.py "$IN" "$ROOT/naive_arms.tsv"
test -f "$ROOT/infercnvpy_arms.tsv" || /usr/bin/time -v $CNV -s $B/run_infercnvpy_bulk.py "$IN" "$ROOT/infercnvpy_arms.tsv" 2> "$ROOT/infercnvpy_time.txt"
$HD -s $B/score.py "$ROOT/score_$(date +%Y%m%d_%H%M%S)" \
  haplodose="$ROOT/haplodose/arms.tsv:median_log_ratio" haplodose_dose="$ROOT/haplodose/arms.tsv:mean_log_dose" \
  naive="$ROOT/naive_arms.tsv" infercnvpy="$ROOT/infercnvpy_arms.tsv"
$HD -s -m pip freeze > "$ROOT/pip_freeze.txt"
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
