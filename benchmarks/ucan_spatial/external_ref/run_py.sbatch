#!/usr/bin/env bash
#SBATCH --job-name=hd_extref_py
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_extref/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_extref/%x_%j.err
# Run one Python script of the external-reference experiment on a compute node.
# Usage: sbatch -J NAME run_py.sbatch SCRIPT [ARGS...]
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s "$@"
