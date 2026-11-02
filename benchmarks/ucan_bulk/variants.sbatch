#!/usr/bin/env bash
#SBATCH --job-name=ucan_bulk_var
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_bulk_var_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_bulk_var_%j.err
# Haplodose bulk variants on U-CAN, scored with the earlier naive and infercnvpy arms.
# Usage: sbatch variants.sbatch INPUT_DIR BASE_ROOT ROOT
set -euo pipefail
IN=$1 BASE=$2 ROOT=$3
mkdir -p "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
run() { local name=$1; shift
  test -f "$ROOT/$name/run.json" || $HD -s -m haplodose.cli bulk --counts "$IN/counts.tsv" --samples "$IN/samples.tsv" \
    --positions "$IN/positions.tsv" --exclude-genes /gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt \
    --out "$ROOT/$name" "$@"; }
run sum_noprog --programmes 0
run median --bin-stat median
run median_noprog --bin-stat median --programmes 0
args=()
for v in sum_noprog median median_noprog; do
  args+=("$v=$ROOT/$v/arms.tsv:median_log_ratio" "${v}_dose=$ROOT/$v/arms.tsv:mean_log_dose")
done
$HD -s benchmarks/ucan_bulk/score.py "$ROOT/score_$(date +%Y%m%d_%H%M%S)" \
  base="$BASE/haplodose/arms.tsv:median_log_ratio" base_dose="$BASE/haplodose/arms.tsv:mean_log_dose" \
  "${args[@]}" naive="$BASE/naive_arms.tsv" infercnvpy="$BASE/infercnvpy_arms.tsv"
$HD -s -m pip freeze > "$ROOT/pip_freeze.txt"
