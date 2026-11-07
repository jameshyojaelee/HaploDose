#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 (STATUS, 2026-11-04 fixtab): rerun ONLY the evaluation rule in an existing gate dir whose evaluation failed
# because fread stops at line 82 of hg38_gencode_v27.txt. The working-copy workflow.sm's evaluation rules are pointed at
# hg38_gencode_v27_fixtab_20261104.txt (method rules unchanged). A dry run must plan exactly one job, the evaluation
# rule; otherwise the job stops before running anything. Method outputs are not touched.
# Usage: sbatch --output=OUT/rerun_eval_%j.out rerun_eval.sbatch OUT DS
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
export RETICULATE_PYTHON=$OUT/pyshim/bin/python
case $DS in
  COLO320) TARGET=results/output_COLO320/evaluation/evaluation_cnv_prediction_auc_cnvkit.txt; RULE=evaluate_cnv_prediction_cnvkit_COLO320; SUF=_cnvkit ;;
  BCC06|BCC06post|MM) TARGET=results/output_$DS/evaluation/evaluation_cnv_prediction_auc_wes.txt; RULE=evaluate_cnv_prediction_wes; SUF=_wes ;;
  MCF7) TARGET=results/output_MCF7/evaluation/evaluation_cnv_prediction_auc.txt; RULE=evaluate_cnv_prediction_wgs; SUF= ;;
  *) echo "unknown dataset $DS"; exit 2 ;;
esac
test ! -e $OUT/pipeline/$TARGET
test ! -e $OUT/gate_comparison.tsv
R=$OUT/rerun_eval_$SLURM_JOB_ID; mkdir -p $R
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ndataset\t$DS\nrule\t$RULE\nstart\t$(date -Is)" > $R/job.tsv
W=$OUT/pipeline/workflow.sm
if ! grep -q hg38_gencode_v27_fixtab_20261104.txt $W; then
  cp $W $W.pre_fixtab
  sed -i 's#infercnv_gene_pos = "data/annotations/hg38_gencode_v27.txt"#infercnv_gene_pos = "data/annotations/hg38_gencode_v27_fixtab_20261104.txt"#' $W
  diff $W.pre_fixtab $W > $OUT/workflow_fixtab.diff || true
fi
test "$(grep -c 'hg38_gencode_v27_fixtab_20261104.txt' $W)" = 14
A=$OUT/pipeline/data/annotations
test -e $A/hg38_gencode_v27_fixtab_20261104.txt || ln -s $H/external/annotations/hg38_gencode_v27_fixtab_20261104.txt $A/hg38_gencode_v27_fixtab_20261104.txt
cd $OUT/pipeline
$MM run -n schmid_r41 snakemake -s workflow.sm --cores 4 -n --rerun-triggers mtime --rerun-incomplete $TARGET > $R/dryrun.log 2>&1
awk '/^Job stats:/{f=1;next} f&&/^$/{f=0} f&&$2~/^[0-9]+$/&&$1!="total"{print $1,$2}' $R/dryrun.log | sort -u > $R/planned_jobs.txt
cat $R/planned_jobs.txt
test "$(cat $R/planned_jobs.txt)" = "$RULE 1" || { echo "DRYRUN_PLANS_MORE_THAN_EVAL; nothing run"; exit 3; }
rc=0
$MM run -n schmid_r41 snakemake -s workflow.sm --cores 4 --resources mem_mb=60000 --rerun-triggers mtime \
  --rerun-incomplete --printshellcmds $TARGET > $R/snakemake.log 2>&1 || rc=$?
if [ $rc -ne 0 ]; then echo -e "snakemake\tnonzero_exit_$rc" >> $R/job.tsv; echo EVAL_FAILED; exit 1; fi
$MM run -n schmid_r41 python -I $S/gate/compare_published.py $OUT/pipeline/results/output_$DS/evaluation $DS "$SUF" \
  $S/gate/published_values.tsv $OUT > $R/compare_published.log 2>&1
$MM env export -n schmid_r41 > $R/env_export.yml
echo -e "end\t$(date -Is)" >> $R/job.tsv
echo RERUN_EVAL_DONE
