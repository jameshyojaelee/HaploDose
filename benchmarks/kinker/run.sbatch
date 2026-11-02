#!/usr/bin/env bash
#SBATCH --job-name=kinker_run
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020 benchmark: Haplodose on both rotations (prepare.py), then evaluate.py vs DepMap.
# Usage: sbatch run.sbatch INPUT_DIR ROOT   (HD_ARGS env: extra haplodose.cli run arguments)
set -euo pipefail
IN=$1 ROOT=$2
mkdir -p "$ROOT"
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
for r in 0 1; do
  test -f "$ROOT/haplodose_r$r/run.json" || /usr/bin/time -v $PY -s -m haplodose.cli run --counts "$IN/counts.h5" \
    --cells "$IN/cells_r$r.tsv" --positions "$IN/positions.tsv" --out "$ROOT/haplodose_r$r" ${HD_ARGS:-} 2> "$ROOT/haplodose_r${r}_time.txt"
done
$PY -s benchmarks/kinker/evaluate.py "$IN" "$ROOT/eval_$(date +%Y%m%d_%H%M%S)" "$ROOT/haplodose_r0" "$ROOT/haplodose_r1"
$PY -s -m pip freeze > "$ROOT/pip_freeze.txt"
