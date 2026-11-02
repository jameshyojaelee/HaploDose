#!/usr/bin/env bash
# Usage: summarise_compare.sh EXTREF_ROOT BASELINE_ROOT (run via run_py.sbatch's environment by sbatch --wrap)
set -euo pipefail
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
HD=/gpfs/commons/home/jameslee/Haplodose/benchmarks/ucan_spatial
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
$PY -s $HD/summarise.py "$1" dna | tee "$1/summarise_dna.txt"
$PY -s $HD/external_ref/compare.py "$2" "$1" | tee "$1/comparison.txt"
$PY -s -m pip freeze > "$1/pip_freeze.txt"
