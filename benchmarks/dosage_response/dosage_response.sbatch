#!/usr/bin/env bash
#SBATCH --job-name=hd_dosage
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/dosage/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/dosage/%x_%j.err
# Per-gene RNA dosage response (U-CAN bulk + FACETS; TCGA COAD/READ replication). Inputs are read only.
# Usage: sbatch dosage_response.sbatch OUTDIR   (OUTDIR must not exist)
set -euo pipefail
OUT=$1
HD=/gpfs/commons/home/jameslee/Haplodose
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir "$OUT"
mkdir "$OUT/scripts"
cp "$HD/benchmarks/dosage_response/dosage_response.py" "$HD/benchmarks/dosage_response/dosage_response.sbatch" "$OUT/scripts/"
$PY -m pip freeze > "$OUT/pip_freeze.txt"
echo "job $SLURM_JOB_ID node $(hostname) $(date -Is)" > "$OUT/job.txt"
$PY -s "$HD/benchmarks/dosage_response/dosage_response.py" "$OUT" --threads "$SLURM_CPUS_PER_TASK"
