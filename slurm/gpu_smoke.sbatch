#!/usr/bin/env bash
#SBATCH --job-name=haplodose_gpu
#SBATCH --account=innovation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/gpu_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/gpu_%j.err
# GPU smoke: same code path with JAX on CUDA (cuda12 plugin in the haplodose env), then the scale
# benchmark on GPU. NOT SUBMITTED: each GPU run needs explicit approval (CRC rule).
# Usage: sbatch slurm/gpu_smoke.sbatch OUT [--cells N] [--specimens S]
set -euo pipefail
export PYTHONNOUSERSITE=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
$PY -s -c "import jax; print(jax.default_backend(), jax.devices())"
cd /gpfs/commons/home/jameslee/Haplodose && $PY -s -m pytest -q tests/test_hmm.py tests/test_likelihood.py
$PY -s benchmarks/scale.py "$@"
