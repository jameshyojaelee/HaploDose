#!/usr/bin/env bash
#SBATCH --job-name=haplodose_wgrid
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/grid/wgrid_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/grid/wgrid_%j.err
# Usage: sbatch slurm/window_grid.sbatch OUT [--quick]
set -euo pipefail
export PYTHONNOUSERSITE=1 JAX_PLATFORMS=cpu
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  /gpfs/commons/home/jameslee/Haplodose/benchmarks/window_grid.py "$@"
