#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 scoring, all lines in one job (team-lead 2026-11-04): BCC06, BCC06post, COLO320, MM (primary) and MCF7
# (line overlap), for HaploDose arm vnext (primary) and sens1424 (descriptive). Per dataset and arm:
#   combined  = their evaluate_results.R patched to add HaploDose (common bins across all methods; primary)
#   permethod = their evaluate_results_permethod.R patched likewise (each method on its own bins; sensitivity)
# No ranking is computed here. Gate dirs and HaploDose dirs are only linked, never written.
# Gate handling per STATUS claim 5 (d)/(g), in gate_check.R: the job stops only if no primary run has every gate row
# within 0.05. A comparator output that misses on any row is still scored; its gap is in gate_check.tsv and every
# scorer table gets a gate_miss column (<table>_gateflag.tsv, flag_tables.R). A run without a gate table, or without
# HaploDose outputs, is not scored and is listed in skipped_runs.tsv; the other runs go on. No overrides.
# HaploDose must have run on byte-identical inputs to the gate; a run that did not is skipped (skipped_runs.tsv).
# INPUTS_TSV columns: run_id (unique, names the output dir), dataset (pipeline name), role, gate_dir, hd_vnext, hd_sens1424.
# Usage: sbatch --account=innovation --qos=innovation -J schmid_score_claim5 --output=OUT/score_%j.out \
#          run_score_all.sbatch OUT [INPUTS_TSV]
set -euo pipefail
H=/gpfs/commons/home/jameslee/Haplodose
# Helpers come from CODE_DIR (copy made at submission by submit_frozen.sh), else the live tree
S=${CODE_DIR:-$H/benchmarks/schmid}
OUT=$1 MAP=${2:-$S/score/score_inputs.tsv}
SRC=$H/external/benchmark_scrnaseq_cnv_callers
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
test "$(git -C $SRC rev-parse HEAD)" = 319ea34ddfc46bdd4c0cb1ed534f4ece0f5c44ca
test ! -e $OUT/vnext && test ! -e $OUT/sens1424
mkdir -p $OUT
cp $MAP $OUT/score_inputs.tsv
echo -e "host\t$(hostname)\njob\t${SLURM_JOB_ID:-none}\nstart\t$(date -Is)" > $OUT/job.tsv
test "$(tail -n +2 $MAP | cut -f1 | sort | uniq -d | wc -l)" = 0 || { echo "DUPLICATE_RUN_ID"; exit 3; }
md5sum $S/evaluate_results_haplodose.R $S/evaluate_results_permethod_haplodose.R $S/lib_haplodose.R \
  $S/snakemake_shim.R $S/score/score_one.R $SRC/snakemake_pipeline/scripts/lib.R > $OUT/scorer_md5.txt

# 1. Gate handling (stops only if no primary run passes every row)
$MM run -n schmid_r41 Rscript $S/score/gate_check.R $MAP $H $OUT > $OUT/gate_check.log 2>&1 || { tail -3 $OUT/gate_check.log; echo "GATE_CHECK_STOP"; exit 3; }
echo -e "run_id\treason" > $OUT/skipped_runs.tsv
awk -F'\t' 'NR>1 && $4=="gate_table_missing"{print $1"\tgate_table_missing"}' $OUT/gate_check.tsv >> $OUT/skipped_runs.tsv

# 2. Preconditions for each run that has a gate table, before anything is scored
: > $OUT/runs_ready.txt
tail -n +2 $MAP | while IFS=$'\t' read -r RID DS ROLE G V SN; do
  grep -qx "$RID" $OUT/runs_to_score.txt || continue
  ok=1
  for A in "$V" "$SN"; do
    for f in ${DS}_haplodose_gene_by_clone.tsv ${DS}_haplodose_cell_clone.tsv; do
      test -s $H/$A/haplodose/$f || { echo -e "$RID\thaplodose_output_missing:$A/$f" >> $OUT/skipped_runs.tsv; ok=0; }
    done
    if [ $ok = 0 ]; then break; fi
    # HaploDose must have run on the same three input files as the gate
    for f in count_matrix.txt sample_annotation.txt ref_groups.txt; do
      x=$(grep "/$f\$" $H/$A/input_md5.txt 2>/dev/null | cut -d' ' -f1 || true); y=$(md5sum < $(readlink -f $H/$G/pipeline/data/input_$DS/$f) | cut -d' ' -f1)
      if [ -z "$x" ] || [ "$x" != "$y" ]; then echo -e "$RID\thaplodose_input_not_identical_to_gate:$A/$f" >> $OUT/skipped_runs.tsv; ok=0; break; fi
    done
    if [ $ok = 0 ]; then break; fi
  done
  if [ $ok = 1 ]; then echo "$RID" >> $OUT/runs_ready.txt; fi
done
test -s $OUT/runs_ready.txt || { echo "NO_RUN_READY"; exit 3; }

# 3. Score, then flag
tail -n +2 $MAP | while IFS=$'\t' read -r RID DS ROLE G V SN; do
  grep -qx "$RID" $OUT/runs_ready.txt || continue
  for ARM in vnext sens1424; do
    [ $ARM = vnext ] && A=$V || A=$SN
    W=$OUT/$ARM/$RID/pipeline; GP=$H/$G/pipeline
    mkdir -p $W/scripts $W/data $W/results/output_$DS $W/logs
    cp $SRC/snakemake_pipeline/scripts/lib.R $W/scripts/
    cp $S/evaluate_results_haplodose.R $S/evaluate_results_permethod_haplodose.R $S/lib_haplodose.R $W/scripts/
    for d in $GP/data/*; do ln -s $d $W/data/$(basename $d); done
    for m in infercnv casper copykat scevan numbat CONICSmat; do ln -s $GP/results/output_$DS/$m $W/results/output_$DS/$m; done
    ln -s $H/$A/haplodose $W/results/output_$DS/haplodose
    echo -e "$RID\t$DS\t$ARM\t$ROLE\t$G\t$A" >> $OUT/score_map.tsv
    for MODE in combined permethod; do
      (cd $W && /usr/bin/time -v $MM run -n schmid_r41 Rscript $S/score/score_one.R $S/snakemake_shim.R $DS $MODE) \
        > $OUT/$ARM/$RID/score_${MODE}.out 2>&1 || { echo "SCORE_FAILED $RID $ARM $MODE"; exit 1; }
    done
    $MM run -n schmid_r41 Rscript $S/score/flag_tables.R $OUT/gate_methods.tsv $RID $W/results/output_$DS/evaluation_haplodose \
      > $OUT/$ARM/$RID/flag_tables.out 2>&1 || { echo "FLAG_FAILED $RID $ARM"; exit 1; }
  done
done
md5sum $S/score/gate_check.R $S/score/flag_tables.R >> $OUT/scorer_md5.txt
echo -e "code_dir\t$S" >> $OUT/job.tsv
$MM env export -n schmid_r41 > $OUT/env_export.yml
echo -e "end\t$(date -Is)" >> $OUT/job.tsv
echo SCORE_ALL_DONE
