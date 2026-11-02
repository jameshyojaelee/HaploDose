#!/usr/bin/env bash
#SBATCH --job-name=kinker_implant
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker subclone design: implant.py on the v2 design, then Haplodose (both rotations) and infercnvpy
# (per pool and rotation, --rule all: every query cell clustered), then score_implant.py.
# Usage: sbatch run_implant.sbatch INPUT_DIR ROOT [GAIN_MODE]   (ROOT/input holds the implanted design)
# HD_ARGS (env): extra haplodose.cli run arguments, e.g. HD_ARGS='--event-windows'.
set -euo pipefail
IN=$1 ROOT=$2 GAIN_MODE=${3:-depth_neutral}
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
mkdir -p "$ROOT"
I=$ROOT/input
test -f "$I/implants.tsv" || $HD -s benchmarks/kinker/implant.py "$IN" "$I" "$GAIN_MODE"
RUNS=()
for r in 0 1; do
  test -f "$ROOT/haplodose_r$r/run.json" || /usr/bin/time -v $HD -s -m haplodose.cli run --counts "$I/counts_r$r.h5" \
    --cells "$I/cells_r$r.tsv" --positions "$I/positions.tsv" --out "$ROOT/haplodose_r$r" ${HD_ARGS:-} 2> "$ROOT/haplodose_r${r}_time.txt"
  RUNS+=("haplodose=$ROOT/haplodose_r$r")
  mkdir -p "$ROOT/cells"
  for p in $(awk -F"\t" 'NR>1 && $11=="single_line" {print $2}' "$I/cells_r$r.tsv" | sort -u); do
    C=$ROOT/cells/r${r}_pool$p.tsv
    test -f "$C" || $HD -s -c "
import pandas as pd
c = pd.read_csv('$I/cells_r$r.tsv', sep='\t', dtype={'pool': str})
o = c.pool.ne('$p')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C', sep='\t', index=False)"
    ls -d "$ROOT/r${r}_pool$p"/infercnvpy_* >/dev/null 2>&1 || { mkdir -p "$ROOT/r${r}_pool$p"; \
      $CNV -s benchmarks/comparators/run_infercnvpy.py --counts "$I/counts_r$r.h5" --cells "$C" \
        --positions "$I/positions.tsv" --out-parent "$ROOT/r${r}_pool$p" --rule all; }
    RUNS+=("infercnvpy=$(ls -d "$ROOT/r${r}_pool$p"/infercnvpy_* | tail -1)")
  done
done
$HD -s benchmarks/kinker/score_implant.py "$I" "$ROOT/score_$(date +%Y%m%d_%H%M%S)" "${RUNS[@]}"
$HD -s -m pip freeze > "$ROOT/pip_freeze_haplodose.txt"
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
