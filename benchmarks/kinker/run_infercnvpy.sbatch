#!/usr/bin/env bash
#SBATCH --job-name=kinker_infercnvpy
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020: infercnvpy 0.6.1 (benchmarks/comparators/run_infercnvpy.py), one run per pool and rotation so
# each pool's queries use only that pool's reference lines (infercnvpy has no per-dataset reference); then
# evaluate.py --tool infercnvpy over all runs.  Usage: sbatch run_infercnvpy.sbatch INPUT_DIR ROOT
set -euo pipefail
IN=$1 ROOT=$2
mkdir -p "$ROOT/cells"
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
for r in 0 1; do
  for p in $(awk -F"\t" 'NR>1 && $11=="single_line" {print $2}' "$IN/cells_r$r.tsv" | sort -u); do
    C=$ROOT/cells/r${r}_pool$p.tsv
    test -f "$C" || $HD -s -c "
import pandas as pd, sys
c = pd.read_csv('$IN/cells_r$r.tsv', sep='\t', dtype={'pool': str})
o = c.pool.ne('$p')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C', sep='\t', index=False)"
    ls -d "$ROOT/r${r}_pool$p"/infercnvpy_* >/dev/null 2>&1 || { mkdir -p "$ROOT/r${r}_pool$p"; \
      $CNV -s benchmarks/comparators/run_infercnvpy.py --counts "$IN/counts.h5" --cells "$C" \
        --positions "$IN/positions.tsv" --out-parent "$ROOT/r${r}_pool$p"; }
  done
done
$HD -s benchmarks/kinker/evaluate.py --tool infercnvpy "$IN" "$ROOT/eval_$(date +%Y%m%d_%H%M%S)" "$ROOT"/r*_pool*/infercnvpy_*
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
