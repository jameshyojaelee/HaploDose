#!/usr/bin/env bash
#SBATCH --job-name=kinker_infercnv_r
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --constraint=v6
# Hard constraint, not --prefer: the R 4.4.1 module died with SIGILL on a v4 Intel node (job 22598463).
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020: R inferCNV (benchmarks/comparators/run_infercnv_r.py), one run per pool and rotation, as
# run_infercnvpy.sbatch; then evaluate.py --tool infercnv_r. Usage: sbatch run_infercnv_r.sbatch INPUT_DIR ROOT
set -euo pipefail
IN=$1 ROOT=$2
mkdir -p "$ROOT/cells"
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
# R env as in benchmarks/comparators/run_one.sbatch.
L=$(env -i HOME=$HOME bash -lc 'module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1; echo "$(command -v Rscript)|$LD_LIBRARY_PATH"')
export RSCRIPT=${L%%|*} R_LD_LIBRARY_PATH=/nfs/sw/easybuild/software/GCCcore/13.2.0/lib64:$(echo "${L#*|}" | tr ':' '\n' | grep -v 'GCCcore/1[12]\.' | paste -sd:):/nfs/sw/easybuild/software/JAGS/4.3.2-foss-2023a/lib
module purge 2>/dev/null || true
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
for r in 0 1; do
  for p in $(awk -F"\t" 'NR>1 && $11=="single_line" {print $2}' "$IN/cells_r$r.tsv" | sort -u); do
    C=$ROOT/cells/r${r}_pool$p.tsv
    test -f "$C" || $HD -s -c "
import pandas as pd
c = pd.read_csv('$IN/cells_r$r.tsv', sep='\t', dtype={'pool': str})
o = c.pool.ne('$p')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C', sep='\t', index=False)"
    ls -d "$ROOT/r${r}_pool$p"/infercnv_r_* >/dev/null 2>&1 || { mkdir -p "$ROOT/r${r}_pool$p"; \
      $CNV -s benchmarks/comparators/run_infercnv_r.py --counts "$IN/counts.h5" --cells "$C" \
        --positions "$IN/positions.tsv" --out-parent "$ROOT/r${r}_pool$p"; }
  done
done
$HD -s benchmarks/kinker/evaluate.py --tool infercnv_r "$IN" "$ROOT/eval_$(date +%Y%m%d_%H%M%S)" "$ROOT"/r*_pool*/infercnv_r_*
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
