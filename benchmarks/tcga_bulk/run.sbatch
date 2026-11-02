#!/usr/bin/env bash
#SBATCH --job-name=tcga_bulk
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/tcga_bulk_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/tcga_bulk_%j.err
# TCGA pan-cancer bulk RNA benchmark against ABSOLUTE: Haplodose bulk (defaults: sum bins, 10 programmes,
# the U-CAN variant winner), naive log ratio and infercnvpy, per cancer type with >= 10 solid tissue normals.
# Plus a raw-count check on COAD+READ (GDC STAR counts vs Xena pseudo-counts, same samples).
# Usage: sbatch run.sbatch ROOT   (ROOT = results/tcga_bulk_<stamp>, must not exist)
set -euo pipefail
ROOT=$1
test ! -e "$ROOT"
mkdir -p "$ROOT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
U=benchmarks/ucan_bulk T=benchmarks/tcga_bulk
MASK=/gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt
IN=$ROOT/input
/usr/bin/time -v $HD -s $T/make_input.py "$IN" --coadread-raw 2> "$ROOT/make_input_time.txt"
tm() { /usr/bin/time -f "%e\t%M" -o "$1" "${@:2}"; }   # wall seconds, peak RSS kB
run_type() { local d=$1
  tm "$d/haplodose_time.tsv" $HD -s -m haplodose.cli bulk --counts "$d/counts.tsv" --samples "$d/samples.tsv" \
    --positions "$d/positions.tsv" --exclude-genes $MASK --out "$d/haplodose" > "$d/haplodose.log" 2>&1
  tm "$d/naive_time.tsv" $HD -s $U/naive.py "$d" "$d/naive_arms.tsv"
  tm "$d/infercnvpy_time.tsv" $CNV -s $U/run_infercnvpy_bulk.py "$d" "$d/infercnvpy_arms.tsv" > "$d/infercnvpy.log" 2>&1
}
types=()
for d in "$IN"/*/; do
  t=$(basename "$d"); run_type "${d%/}"
  [[ $t == coadread_* ]] || types+=("$t=${d%/}")
done
M=(haplodose=haplodose/arms.tsv:median_log_ratio haplodose_dose=haplodose/arms.tsv:mean_log_dose
   naive=naive_arms.tsv infercnvpy=infercnvpy_arms.tsv)
$HD -s $T/score.py "$ROOT/score" "${types[@]}" --methods "${M[@]}"
$HD -s $T/score.py "$ROOT/score_coadread_counts" xena=$IN/coadread_xena raw=$IN/coadread_raw --methods "${M[@]}"
$HD -s -m pip freeze > "$ROOT/pip_freeze.txt"
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
