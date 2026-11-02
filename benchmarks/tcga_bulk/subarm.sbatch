#!/usr/bin/env bash
#SBATCH --job-name=tcga_subarm
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/tcga_subarm_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/tcga_subarm_%j.err
# Sub-arm (10 Mb gene-anchored windows) TCGA benchmark on an existing tcga_bulk run (no Haplodose re-profiling).
# Usage: sbatch subarm.sbatch RUN_ROOT OUT   (OUT = RUN_ROOT/subarm_<stamp>, must not exist)
set -euo pipefail
RUN=$1 OUT=$2
test ! -e "$OUT"
mkdir -p "$OUT/icnv"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
types=()
for d in "$RUN"/input/*/; do d=${d%/}; t=$(basename "$d")
  [[ $t == coadread_* ]] && continue
  $CNV -s benchmarks/tcga_bulk/subarm_infercnvpy.py "$d" "$OUT/icnv/$t" > "$OUT/icnv/$t.log" 2>&1
  types+=("$t=$d:$OUT/icnv/$t")
done
/usr/bin/time -v $HD -s benchmarks/tcga_bulk/subarm.py "$OUT/score" "${types[@]}" 2> >(grep -v -i warn > "$OUT/subarm_stderr.txt")
$HD -s -m pip freeze > "$OUT/pip_freeze.txt"
$CNV -s -m pip freeze > "$OUT/pip_freeze_scrna_cnv.txt"
