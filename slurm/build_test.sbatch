#!/usr/bin/env bash
#SBATCH --job-name=haplodose_build_test
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/build_test_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/build_test_%j.err
set -euo pipefail
export FREEZE_OUT=logs/pip_freeze_${SLURM_JOB_ID}.txt
scripts/build_and_test.sh "$@"
