#!/usr/bin/env bash
#SBATCH --job-name=hd_ycc
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ycc/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ycc/%x_%j.err
# Haplodose with alleles on one GSE273567 (YCC) sample, CRC inputs read only; scored against native Numbat
# (no DNA truth). YCC097/155/162 have no positions.tsv: YCC230's is used (same Cell Ranger GRCh38-2020-A).
# Usage: sbatch -J hd_ycc_S run.sbatch SAMPLE ROOT
set -euo pipefail
S=$1 ROOT=$2
C=/gpfs/commons/home/jameslee/CRC
HD=/gpfs/commons/home/jameslee/Haplodose
if [[ -d $C/results/cnv_caller_benchmark/$S/input ]]; then IN=$C/results/cnv_caller_benchmark/$S/input; NB=$C/results/numbat_11p15/$S
else IN=$C/results/numbat_gse273567_extra/$S/caller_input; NB=$C/results/numbat_gse273567_extra/$S; fi
POS=$IN/positions.tsv; [[ -f $POS ]] || POS=$C/results/cnv_caller_benchmark/YCC230/input/positions.tsv
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
[[ -f $ROOT/$S/haplodose/qc.json ]] ||
  $PY -s -m haplodose.cli run --counts $C/data/external/gse273567/$S/alignment/$S/outs/filtered_feature_bc_matrix.h5 \
      --cells $IN/cells.tsv --positions $POS --exclude-genes $C/src/44_haplodose/masks.txt \
      --alleles $NB/alleles/${S}_allele_counts.tsv.gz --allele-exclude-genes $C/src/44_haplodose/allele_exclude_genes.txt \
      --allele-exclude-sites $C/src/44_haplodose/allele_exclude_sites.txt --out $ROOT/$S/haplodose
$PY -s $HD/benchmarks/comparators/ycc_vs_numbat.py $ROOT/$S/haplodose $NB/numbat $S $ROOT/$S/vs_numbat
$PY -s -m pip freeze > $ROOT/$S/pip_freeze.txt
