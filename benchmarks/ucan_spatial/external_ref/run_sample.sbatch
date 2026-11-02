#!/usr/bin/env bash
#SBATCH --job-name=hd_extref
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_extref/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_extref/%x_%j.err
# One U-CAN Visium HD section with the external normal-epithelium reference added: build the combined
# input (section input from CRC, read only, plus ROOT/reference), run Haplodose, score against WGS FACETS.
# Copy of ../run_sample.sbatch without the CRC prepare step.
# Usage: sbatch -J hd_extref_U0007 run_sample.sbatch SAMPLE ROOT [extra haplodose run flags]
set -euo pipefail
S=$1 ROOT=$2
shift 2
HD=/gpfs/commons/home/jameslee/Haplodose
CRC=/gpfs/commons/home/jameslee/CRC
OUT=$ROOT/$S
mkdir -p "$OUT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=$CRC:/gpfs/commons/groups:/gpfs/commons/datasets
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
[[ -f $OUT/input/input.json ]] || { rm -rf "$OUT/input.tmp"; $PY -s $HD/benchmarks/ucan_spatial/external_ref/make_input.py \
    "$CRC/results/cnv_caller_benchmark/$S/input" "$ROOT/reference" "$OUT/input.tmp" && mv -T "$OUT/input.tmp" "$OUT/input"; }
IN=$OUT/input
[[ -f $OUT/haplodose/qc.json ]] ||
  $PY -s -m haplodose.cli run --counts "$IN/counts.h5" --cells "$IN/cells.tsv" --positions "$IN/positions.tsv" \
      --exclude-genes $CRC/src/44_haplodose/masks.txt --out "$OUT/haplodose" "$@"
if [[ ! -f $OUT/dna/summary.json ]]; then
  rmdir "$OUT/dna" 2>/dev/null || true
  $PY -s $HD/benchmarks/ucan_spatial/dna_truth.py "$OUT/haplodose" "$S" "$OUT/dna"
fi
$PY -s -m pip freeze > "$OUT/pip_freeze.txt"
