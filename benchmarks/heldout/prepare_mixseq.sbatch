#!/usr/bin/env bash
#SBATCH --job-name=prepare_mixseq
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout/%x_%j.out
# Usage: sbatch prepare_mixseq.sbatch OUT
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s -I benchmarks/heldout/prepare_mixseq.py \
  data/external/mixseq_mcfarland2020/extracted "$1"
