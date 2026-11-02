#!/usr/bin/env bash
#SBATCH --job-name=scaling_prepare
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/scaling/%x_%j.out
# Usage: sbatch prepare.sbatch OUT_DIR   (input sizes for the compute scaling benchmark)
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s benchmarks/scaling/prepare.py \
  results/kinker_input_v2_20261101_1551 "$1/input" 5000 20000 50000 100000 1000000
