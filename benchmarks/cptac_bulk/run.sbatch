#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# CPTAC-3 bulk sub-arm benchmark (prespecified claim 6, docs/STATUS.md). One stage per job:
#   input        make_input.py DATA -> ROOT/input (CPTAC Solid Tissue Normal panel per cohort)
#   group GROUP  frozen 14:24 HaploDose bulk (same flags as the TCGA run) + infercnvpy windows for that group
#   score        score.py (subarm.py with hg38 AscatNGS truth + the primary test) -> ROOT/score. Only after the lead confirms.
# HaploDose code: results/heldout_20261102_1424/code (read-only), hash-checked at start of every stage.
# Evaluation code hashes are written to ROOT/eval_code_sha256.txt by the submit step and checked here.
# Usage: sbatch -J cptac_<stage> --output=ROOT/logs/%x_%j.out run.sbatch ROOT STAGE [GROUP]
# DATA (md5-checked GDC download) is fixed below.
set -euo pipefail
ROOT=$1 STAGE=$2 GROUP=${3:-}
F=/gpfs/commons/home/jameslee/Haplodose/results/heldout_20261102_1424
(cd $F && sha256sum --quiet -c code_sha256.txt)
sha256sum --quiet -c "$ROOT/eval_code_sha256.txt"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
MASK=/gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt
DATA=/gpfs/commons/home/jameslee/Haplodose/data/external/cptac3_20261103_0616
case $STAGE in
input)
  /usr/bin/time -v $HD -s benchmarks/cptac_bulk/make_input.py "$DATA" "$ROOT/input" 2> "$ROOT/logs/make_input_time.txt"
  $HD -s -m pip freeze > "$ROOT/pip_freeze_haplodose.txt"
  $CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt" ;;
group)
  d=$ROOT/input/$GROUP
  test -f "$d/haplodose/samples.tsv" || PYTHONPATH=$F/code /usr/bin/time -f "%e\t%M" -o "$d/haplodose_time.tsv" \
    $HD -s -m haplodose.cli bulk --counts "$d/counts.tsv" --samples "$d/samples.tsv" --positions "$d/positions.tsv" \
    --exclude-genes $MASK --out "$d/haplodose" > "$d/haplodose.log" 2>&1
  mkdir -p "$ROOT/icnv"
  test -f "$ROOT/icnv/$GROUP/infercnvpy_windows.npz" || $CNV -s benchmarks/tcga_bulk/subarm_infercnvpy.py "$d" \
    "$ROOT/icnv/$GROUP" > "$ROOT/icnv/$GROUP.log" 2>&1 ;;
score)
  /usr/bin/time -v $HD -s benchmarks/cptac_bulk/score.py "$ROOT/score" "$ROOT/input" "$ROOT/icnv" \
    2> >(grep -v -i warn > "$ROOT/logs/score_stderr.txt") ;;
*) echo "unknown stage $STAGE"; exit 2 ;;
esac
