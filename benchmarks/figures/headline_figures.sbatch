#!/usr/bin/env bash
#SBATCH --job-name=hd_headline_figs
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Usage: sbatch -o <new_output_dir>/slurm_%j.out benchmarks/figures/headline_figures.sbatch <new_output_dir>
set -euo pipefail
/gpfs/commons/home/jameslee/micromamba/envs/cnv_field/bin/Rscript benchmarks/figures/headline_figures.R "$1"
