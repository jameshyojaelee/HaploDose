#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 (STATUS 2026-11-04 ~16:45, MCF7 Numbat failure): evaluation only, in an existing gate dir whose run_numbat
# failed. Steps: point the 14 evaluation gene_pos inputs at the fixtab file (if not done), append a copy of RULE with only
# numbat_obj removed (make_nonumbat_rule.py; output evaluation_nonumbat/), dry run must plan exactly "<RULE>_nonumbat 1"
# (else exit 3, nothing run), run it, then compare_published.py so Numbat rows appear with ours = NA.
# Usage: sbatch --output=OUT/rerun_eval_nonumbat_%j.out rerun_eval_nonumbat.sbatch OUT DS   (DS = MCF7 only for now)
set -euo pipefail
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
export RETICULATE_PYTHON=$OUT/pyshim/bin/python
case $DS in
  MCF7) RULE=evaluate_cnv_prediction_wgs; SUF= ;;
  *) echo "no-Numbat evaluation is approved for MCF7 only"; exit 2 ;;
esac
TARGET=results/output_$DS/evaluation_nonumbat/evaluation_cnv_prediction_auc$SUF.txt
test ! -e $OUT/pipeline/$TARGET
test ! -e $OUT/gate_comparison.tsv
R=$OUT/rerun_eval_nonumbat_$SLURM_JOB_ID; mkdir -p $R
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ndataset\t$DS\nrule\t${RULE}_nonumbat\ncode_dir\t$S\nstart\t$(date -Is)" > $R/job.tsv
W=$OUT/pipeline/workflow.sm
if ! grep -q hg38_gencode_v27_fixtab_20261104.txt $W; then
  cp $W $W.pre_fixtab
  sed -i 's#infercnv_gene_pos = "data/annotations/hg38_gencode_v27.txt"#infercnv_gene_pos = "data/annotations/hg38_gencode_v27_fixtab_20261104.txt"#' $W
  diff $W.pre_fixtab $W > $OUT/workflow_fixtab.diff || true
fi
test "$(grep -c 'hg38_gencode_v27_fixtab_20261104.txt' $W)" = 14
A=$OUT/pipeline/data/annotations
test -e $A/hg38_gencode_v27_fixtab_20261104.txt || ln -s $H/external/annotations/hg38_gencode_v27_fixtab_20261104.txt $A/hg38_gencode_v27_fixtab_20261104.txt
grep -q "^rule ${RULE}_nonumbat:" $W && { echo "nonumbat rule already present; stop"; exit 2; }
cp $W $W.pre_nonumbat
$MM run -n schmid_r41 python -I $S/gate/make_nonumbat_rule.py $W $RULE > $R/nonumbat_rule.txt
diff $W.pre_nonumbat $W > $OUT/workflow_nonumbat.diff || true
cd $OUT/pipeline
$MM run -n schmid_r41 snakemake -s workflow.sm --cores 4 -n --rerun-triggers mtime --rerun-incomplete $TARGET > $R/dryrun.log 2>&1
awk '/^Job stats:/{f=1;next} f&&/^$/{f=0} f&&$2~/^[0-9]+$/&&$1!="total"{print $1,$2}' $R/dryrun.log | sort -u > $R/planned_jobs.txt
cat $R/planned_jobs.txt
test "$(cat $R/planned_jobs.txt)" = "${RULE}_nonumbat 1" || { echo "DRYRUN_PLANS_MORE_THAN_EVAL; nothing run"; exit 3; }
rc=0
$MM run -n schmid_r41 snakemake -s workflow.sm --cores 4 --resources mem_mb=60000 --rerun-triggers mtime \
  --rerun-incomplete --printshellcmds $TARGET > $R/snakemake.log 2>&1 || rc=$?
if [ $rc -ne 0 ]; then echo -e "snakemake\tnonzero_exit_$rc" >> $R/job.tsv; echo EVAL_FAILED; exit 1; fi
$MM run -n schmid_r41 python -I $S/gate/compare_published.py $OUT/pipeline/results/output_$DS/evaluation_nonumbat $DS "$SUF" \
  $S/gate/published_values.tsv $OUT > $R/compare_published.log 2>&1
echo -e "evaluation_dir\tresults/output_$DS/evaluation_nonumbat" >> $R/job.tsv
$MM env export -n schmid_r41 > $R/env_export.yml
echo -e "end\t$(date -Is)" >> $R/job.tsv
echo RERUN_EVAL_NONUMBAT_DONE
