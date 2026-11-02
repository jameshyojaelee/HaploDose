#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --constraint=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout/%x_%j.err
# Held-out evaluation (docs/heldout_protocol.md), one tool on one design. Hard --constraint=v6: results move
# slightly with CPU generation, so every held-out run uses one generation (R inferCNV also needs it).
# Haplodose runs from the frozen copy ROOT/code (its SHA-256 list is checked at start); benchmark and comparator
# scripts run from the repo only after their SHA-256 matches the frozen copy in ROOT/eval_code (else the job stops).
# Amendment 4 (2026-11-02, before any method run): stderr goes to time_<tool>.txt so the run-directory glob no
# longer matches it; only runs with qc.json (the last file every wrapper writes) are scored; a tool that fails on
# a pool, or fails outright, is listed in <tool>_failures.txt and the job exits 0, so the afterany analysis runs
# and scores the failure as the protocol says (E1 NaN, E2 ari_all 0); per-pool cell files are written atomically.
# Usage: sbatch -J ho_<tool>_<design> run_tool.sbatch ROOT TOOL DESIGN   (DESIGN: single | implant)
# HO_INPUT overrides the single-design input dir, only for the pre-run check on Kinker data; unset for held-out runs.
set -euo pipefail
ROOT=$1 TOOL=$2 DESIGN=$3
INPUT=${HO_INPUT:-results/heldout_mixseq_input_20261102_1159}
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
check() {  # repo file, frozen copy
  [[ $(sha256sum < "$1") == $(sha256sum < "$2") ]] || { echo "hash mismatch: $1 vs $2" >&2; exit 2; }
}
for f in evaluate compare_tools score_implant implant; do check benchmarks/kinker/$f.py $ROOT/eval_code/$f.py; done
for f in $ROOT/eval_code/comparators/*; do check benchmarks/comparators/$(basename $f) $f; done
(cd $ROOT && sha256sum --quiet -c code_sha256.txt) || { echo "frozen code hash mismatch in $ROOT" >&2; exit 2; }
if [[ $DESIGN == single ]]; then IN=$INPUT RULE=cluster; else IN=$ROOT/implant/input RULE=all; fi
OUT=$ROOT/$DESIGN
FAIL=$OUT/${TOOL}_failures.txt
mkdir -p $OUT/cells
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ntool\t$TOOL\ndesign\t$DESIGN\nstart\t$(date -Is)" > $OUT/${TOOL}_job.tsv
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
counts() { [[ $DESIGN == single ]] && echo "$IN/counts.h5" || echo "$IN/counts_r$1.h5"; }
RUNS=()
if [[ $TOOL == haplodose ]]; then
  export PYTHONPATH=$ROOT/code JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
  export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
  $HD -s -c "import haplodose; print(haplodose.__file__)" > $OUT/haplodose_code_path.txt
  for r in 0 1; do
    test -f $OUT/haplodose_r$r/qc.json || /usr/bin/time -v $HD -s -m haplodose.cli run --counts "$(counts $r)" \
      --cells "$IN/cells_r$r.tsv" --positions "$IN/positions.tsv" --out "$OUT/haplodose_r$r" 2> "$OUT/haplodose_r${r}_time.txt" \
      || { echo "failed haplodose r$r" | tee -a $FAIL; tail -20 "$OUT/haplodose_r${r}_time.txt"; }
    if [[ -f $OUT/haplodose_r$r/qc.json ]]; then RUNS+=("$OUT/haplodose_r$r"); else echo "incomplete $OUT/haplodose_r$r" >> $FAIL; fi
  done
  unset PYTHONPATH
else
  if [[ $TOOL == infercnv_r ]]; then
    L=$(env -i HOME=$HOME bash -lc 'module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1; echo "$(command -v Rscript)|$LD_LIBRARY_PATH"')
    export RSCRIPT=${L%%|*} R_LD_LIBRARY_PATH=/nfs/sw/easybuild/software/GCCcore/13.2.0/lib64:$(echo "${L#*|}" | tr ':' '\n' | grep -v 'GCCcore/1[12]\.' | paste -sd:):/nfs/sw/easybuild/software/JAGS/4.3.2-foss-2023a/lib
  fi
  unset PYTHONPATH PYTHONHOME R_LIBS R_LIBS_USER
  for r in 0 1; do
    for p in $(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["design"]=="single_line" {print $h["pool"]}' "$IN/cells_r$r.tsv" | sort -u); do
      C=$OUT/cells/r${r}_pool$p.tsv
      # Jobs for the same design share $C: write to a job-private temp file, then rename (atomic on one filesystem).
      test -f "$C" || { $HD -s -c "
import pandas as pd
c = pd.read_csv('$IN/cells_r$r.tsv', sep='\t', dtype={'pool': str})
o = c.pool.ne('$p')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C.tmp.$SLURM_JOB_ID', sep='\t', index=False)" && mv -f "$C.tmp.$SLURM_JOB_ID" "$C"; }
      D=$OUT/r${r}_pool$p
      mkdir -p $D
      ls $D/${TOOL}_2*/qc.json >/dev/null 2>&1 && continue
      ARGS=(--counts "$(counts $r)" --cells "$C" --positions "$IN/positions.tsv" --out-parent "$D" --rule $RULE)
      case $TOOL in
        infercnvpy) $CNV -s benchmarks/comparators/run_infercnvpy.py "${ARGS[@]}" 2> $D/time_${TOOL}.txt ;;
        infercnv_r) $CNV -s benchmarks/comparators/run_infercnv_r.py "${ARGS[@]}" 2> $D/time_${TOOL}.txt ;;
        copykat|scevan) $CNV -s benchmarks/comparators/run_r_field.py --tool $TOOL "${ARGS[@]}" 2> $D/time_${TOOL}.txt ;;
        *) echo "unknown tool $TOOL" >&2; exit 2 ;;
      esac || { echo "failed $TOOL r$r pool$p" | tee -a $FAIL; tail -20 $D/time_${TOOL}.txt; }
    done
  done
  # Completed runs only: run directories are <tool>_<YYYYmmdd_HHMMSS>; a failed pool leaves one without qc.json.
  for d in $OUT/r*_pool*/${TOOL}_2*/; do
    [[ -d $d ]] || continue
    if [[ -f $d/qc.json ]]; then RUNS+=("${d%/}"); else echo "incomplete ${d%/}" >> $FAIL; fi
  done
fi
if [[ $DESIGN == single ]]; then
  if (( ${#RUNS[@]} )); then
    $HD -s benchmarks/kinker/evaluate.py --tool $TOOL "$IN" "$OUT/eval_$TOOL" "${RUNS[@]}" \
      || echo "evaluate failed for $TOOL" | tee -a $FAIL
  else
    echo "no completed runs for $TOOL" | tee -a $FAIL
  fi
fi
echo -e "end\t$(date -Is)" >> $OUT/${TOOL}_job.tsv
