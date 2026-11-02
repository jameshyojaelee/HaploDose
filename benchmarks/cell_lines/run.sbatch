#!/usr/bin/env bash
#SBATCH --job-name=hd_celllines
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/celllines/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/celllines/%x_%j.err
# GSE149224 cell-line benchmark: build input, run Haplodose, score against DepMap. Usage: run.sbatch ROOT
set -euo pipefail
ROOT=$1
HD=/gpfs/commons/home/jameslee/Haplodose
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
[[ -f $ROOT/input/counts.h5 ]] || $PY -s $HD/benchmarks/cell_lines/prepare.py "$ROOT/input"
$PY -s -m haplodose.cli run --counts "$ROOT/input/counts.h5" --cells "$ROOT/input/cells.tsv" \
    --positions "$ROOT/input/positions.tsv" --out "$ROOT/haplodose"
$PY -s $HD/benchmarks/cell_lines/evaluate.py "$ROOT/haplodose" "$ROOT/input" "$ROOT/eval"
$PY -s $HD/benchmarks/cell_lines/bulk_ceiling.py "$ROOT/haplodose" "$ROOT/eval/bulk_ceiling.json"
$PY -s -m pip freeze > "$ROOT/pip_freeze.txt"
