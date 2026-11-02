#!/usr/bin/env bash
#SBATCH --job-name=hd_ucan
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan/%x_%j.err
# One U-CAN Visium HD section: prepare input (CRC adapter, read only) if none exists, run Haplodose
# (spatial mode), compare arms with the patient's WGS FACETS calls.
# Usage: sbatch -J hd_ucan_U0007 run_sample.sbatch SAMPLE ROOT
set -euo pipefail
S=$1 ROOT=$2
HD=/gpfs/commons/home/jameslee/Haplodose
CRC=/gpfs/commons/home/jameslee/CRC
OUT=$ROOT/$S
mkdir -p "$OUT"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
IN=$CRC/results/cnv_caller_benchmark/$S/input
if [[ ! -f $IN/counts.h5 ]]; then
  # CRC's write guard only allows outputs inside the CRC project: stage there, then move into the
  # (new) cnv_caller_benchmark/$S directory, the same layout as the four sections prepared earlier.
  [[ ! -e $CRC/results/cnv_caller_benchmark/$S ]] || { echo "partial $CRC/results/cnv_caller_benchmark/$S exists" >&2; exit 1; }
  STAGE=$CRC/results/.cnv_caller_benchmark.$S.input.stage.$SLURM_JOB_ID
  ( cd "$CRC" && set +u && source config/env.sh && set -u && unset PYTHONPATH PYTHONHOME LD_LIBRARY_PATH &&
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    /gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python -s -B \
      src/43_cnv_caller/spatial.py prepare --sample "$S" --out "$STAGE/input" --gtf "$GTF" --contig-map "$CONTIG_MAP" &&
    record_provenance "$STAGE" haplodose_ucan_input )
  mv -T "$STAGE" "$CRC/results/cnv_caller_benchmark/$S"
fi
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
[[ -f $OUT/haplodose/qc.json ]] ||
  $PY -s -m haplodose.cli run --counts "$IN/counts.h5" --cells "$IN/cells.tsv" --positions "$IN/positions.tsv" \
      --exclude-genes "${EXCLUDE:-$CRC/src/44_haplodose/masks.txt}" --out "$OUT/haplodose"
if [[ ! -f $OUT/dna/summary.json ]]; then
  rmdir "$OUT/dna" 2>/dev/null || true      # empty directory left by a failed comparison
  $PY -s $HD/benchmarks/ucan_spatial/dna_truth.py "$OUT/haplodose" "$S" "$OUT/dna"
fi
$PY -s -m pip freeze > "$OUT/pip_freeze.txt"
