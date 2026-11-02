#!/usr/bin/env bash
# Submit the head-to-head comparator runs, one sbatch job per tool x dataset (no arrays).
# Usage: submit_all.sh ROOT   (ROOT is a new or existing results/comparators_<stamp> directory)
# Skips a tool x dataset whose ROOT/<dataset>/<tool>_* already has qc.json.
set -euo pipefail
ROOT=$1
cd /gpfs/commons/home/jameslee/Haplodose
B=/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark
M=/gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt
IMPRINTED=$ROOT/exclude_imprinted.txt
mkdir -p "$ROOT" logs/comparators
# Imprinted-gene exclusion from ~/CRC/src/42_meeting_followups/T2_run.R (Numbat runs).
printf '%s\n' IGF2 H19 KCNQ1OT1 CDKN1C KCNQ1 PHLDA2 SLC22A18 MEG3 DLK1 MEST PEG3 PEG10 SNRPN SNURF NNAT PLAGL1 GNAS > "$IMPRINTED"
declare -A MEM=([infercnvpy]=16G [prototype]=16G [infercnv_r]=64G [numbat]=96G)

submit() {  # tool dataset scorer runner-args...
  local tool=$1 ds=$2 scorer=$3; shift 3
  if ls "$ROOT/$ds/${tool}_"*/qc.json >/dev/null 2>&1; then echo "skip $tool $ds"; return; fi
  local ex=(--exclude-genes "$M"); [[ $tool == prototype ]] && ex=()   # the prototype masks IGF2 itself
  [[ $tool == numbat ]] && ex=(--exclude-genes "$IMPRINTED")
  sbatch --parsable -J "cmp_${tool}_${ds}" --mem="${MEM[$tool]}" benchmarks/comparators/run_one.sbatch \
    "$tool" "$ROOT/$ds" "$scorer" -- "$@" "${ex[@]}" | sed "s/^/$tool $ds /"
}
contract() { echo --counts "$1/counts.h5" --cells "$1/cells.tsv" --positions "$1/positions.tsv"; }

# (a) U-CAN Visium HD (probe-based; Numbat does not apply).
for IN in "$B"/U*/input; do
  S=$(basename "$(dirname "$IN")")
  [[ -f $IN/counts.h5 ]] || { echo "no input $S"; continue; }
  for t in infercnvpy prototype infercnv_r; do submit $t "$S" "ucan:$S" $(contract "$IN"); done
done
# (b) GSE149224 cell lines.
for t in infercnvpy prototype infercnv_r; do submit $t celllines cell_lines $(contract results/celllines_20261101_0853/input); done
# (c) Moorman scRNA vs MSK-IMPACT.
for t in infercnvpy prototype infercnv_r; do submit $t Moorman moorman $(contract "$B/Moorman_full/input"); done
# (d) YCC with alleles: every tool, Numbat included. Counts are the CellRanger filtered H5.
declare -A ALLELES=([YCC061]=/gpfs/commons/home/jameslee/CRC/results/numbat_11p15/YCC061/alleles/YCC061_allele_counts.tsv.gz
                    [YCC230]=/gpfs/commons/home/jameslee/CRC/results/numbat_11p15/YCC230/alleles/YCC230_allele_counts.tsv.gz)
for S in YCC061 YCC230; do
  A=(--counts /gpfs/commons/home/jameslee/CRC/data/external/gse273567/$S/alignment/$S/outs/filtered_feature_bc_matrix.h5
     --cells "$B/$S/input/cells.tsv" --positions "$B/$S/input/positions.tsv")
  for t in infercnvpy prototype infercnv_r; do submit $t "$S" none "${A[@]}"; done
  submit numbat "$S" none "${A[@]}" --alleles "${ALLELES[$S]}" --specimen "$S"
done
