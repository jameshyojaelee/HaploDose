#!/usr/bin/env bash
#SBATCH --job-name=kinker_field
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020: CopyKAT or SCEVAN (benchmarks/comparators/run_r_field.py, env cnv_field), one run per pool and
# rotation as run_infercnvpy.sbatch, then evaluate.py --tool TOOL (single-line design) or nothing (implant).
# Usage: sbatch run_field.sbatch TOOL INPUT_DIR ROOT [RULE] [POOLS]   (COUNTS env: counts file per rotation,
#   '{r}' replaced by the rotation; default INPUT_DIR/counts.h5. POOLS: space-separated subset, for smoke tests.)
set -euo pipefail
TOOL=$1 IN=$2 ROOT=$3 RULE=${4:-cluster} POOLS=${5:-}
mkdir -p "$ROOT/cells"
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
unset PYTHONPATH PYTHONHOME R_LIBS R_LIBS_USER
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
for r in 0 1; do
  CELLS=$IN/cells_r$r.tsv
  C0=${COUNTS:-$IN/counts.h5}; C0=${C0//\{r\}/$r}
  for p in $(awk -F"\t" 'NR>1 && $11=="single_line" {print $2}' "$CELLS" | sort -u); do
    [[ -z $POOLS || " $POOLS " == *" $p "* ]] || continue
    C=$ROOT/cells/r${r}_pool$p.tsv
    test -f "$C" || $HD -s -c "
import pandas as pd
c = pd.read_csv('$CELLS', sep='\t', dtype={'pool': str})
o = c.pool.ne('$p')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C', sep='\t', index=False)"
    ls -d "$ROOT/r${r}_pool$p"/${TOOL}_* >/dev/null 2>&1 || { mkdir -p "$ROOT/r${r}_pool$p"; \
      /usr/bin/time -v $CNV -s benchmarks/comparators/run_r_field.py --tool $TOOL --counts "$C0" --cells "$C" \
        --positions "$IN/positions.tsv" --out-parent "$ROOT/r${r}_pool$p" --rule $RULE 2> "$ROOT/r${r}_pool${p}_time.txt" \
        || { tail -30 "$ROOT/r${r}_pool${p}_time.txt"; exit 1; }; }
  done
done
if [[ -z $POOLS && -e $IN/depmap_cn.tsv.gz ]]; then
  $HD -s benchmarks/kinker/evaluate.py --tool $TOOL "$IN" "$ROOT/eval_$(date +%Y%m%d_%H%M%S)" "$ROOT"/r*_pool*/${TOOL}_*
fi
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
