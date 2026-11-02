#!/usr/bin/env bash
#SBATCH --job-name=ucan_allele_score
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/%x_%j.out
# Arm allelic imbalance per sample (arm_baf.py), then scoring vs FACETS, WGD and IGF2 TD (score_alleles.py).
# Usage: sbatch score.sbatch PILEUP_DIR SITES_DIR BULK_ARMS_TSV ROOT
set -euo pipefail
PILE=$1 SITES=$2 BULK=$3 ROOT=$4
mkdir -p "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
A=benchmarks/ucan_bulk/alleles
test -f "$ROOT/arms_baf.tsv" || /usr/bin/time -v $PY -s $A/arm_baf.py "$PILE" "$SITES" "$ROOT/arms_baf.tsv" > "$ROOT/arm_baf.log" 2>&1
$PY -s $A/score_alleles.py "$ROOT/arms_baf.tsv" "$BULK" "$ROOT/score_$(date +%Y%m%d_%H%M%S)"
