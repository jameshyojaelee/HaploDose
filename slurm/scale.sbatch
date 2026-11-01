#!/usr/bin/env bash
#SBATCH --job-name=haplodose_scale
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/scale_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/scale_%j.err
# CPU scale benchmark. Usage: sbatch slurm/scale.sbatch OUT [--cells N] [--specimens S]
set -euo pipefail
export PYTHONNOUSERSITE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s /gpfs/commons/home/jameslee/Haplodose/benchmarks/scale.py "$@"
