#!/usr/bin/env bash
#SBATCH --job-name=hd_moorman
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/moorman/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/moorman/%x_%j.err
# Moorman scRNA (CRC prepared input, read only): Haplodose run, then arm comparison with MSK-IMPACT.
# Usage: run.sbatch ROOT
set -euo pipefail
ROOT=$1
IN=/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark/Moorman_full/input
HD=/gpfs/commons/home/jameslee/Haplodose
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
[[ -f $ROOT/haplodose/qc.json ]] ||
  $PY -s -m haplodose.cli run --counts $IN/counts.h5 --cells $IN/cells.tsv --positions $IN/positions.tsv \
      --exclude-genes /gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt --out "$ROOT/haplodose"
$PY -s $HD/benchmarks/moorman/dna_truth.py "$ROOT/haplodose" "$ROOT/dna"
$PY -s -m pip freeze > "$ROOT/pip_freeze.txt"
