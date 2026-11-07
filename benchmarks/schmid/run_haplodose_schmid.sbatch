#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -J masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 HaploDose run on one Schmid dataset (docs/STATUS.md "Claim 5 frozen code, 2026-11-03 ~20:40").
#   ARM=vnext   (primary): results/vnext_freeze_20261103_0744/code/python, flags --split-mode evenodd
#               --split-null-draws 19 --split-null-model cov, --seed 20261103 (CLI) and wrapper split seed 20261103.
#   ARM=sens1424 (descriptive sensitivity, not part of the claim): results/heldout_20261102_1424/code, CLI defaults.
# Both: value log_ratio_gene, reference holdout 0.3, CPU (PLAN.md item 4). Code SHA-256 checked before running.
# No scoring against truth here.
# Usage: sbatch --output=OUT/DS/job_%j.out run_haplodose_schmid.sbatch OUT DS INPUT_DIR ARM
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2 INP=$3 ARM=$4
H=/gpfs/commons/home/jameslee/Haplodose
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
case $ARM in
  vnext)    F=$H/results/vnext_freeze_20261103_0744; (cd $F/code && sha256sum --quiet -c ../code_sha256.txt)
            CODE=$F/code/python; HD="--split-mode evenodd --split-null-draws 19 --split-null-model cov --seed 20261103" ;;
  sens1424) F=$H/results/heldout_20261102_1424; (cd $F && sha256sum --quiet -c code_sha256.txt)
            CODE=$F/code; HD="" ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
D=$OUT/$DS; mkdir -p $D
test ! -e $D/haplodose
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\narm\t$ARM\ncode\t$CODE\nhd_args\t$HD\ncpu\t$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)\nstart\t$(date -Is)" > $D/job.tsv
md5sum $INP/count_matrix.txt $INP/sample_annotation.txt $INP/ref_groups.txt > $D/input_md5.txt
/usr/bin/time -v $PY -I $S/run_haplodose.py \
  --matrix $INP/count_matrix.txt --annot $INP/sample_annotation.txt --ref-groups $INP/ref_groups.txt \
  --gene-pos $H/external/annotations/hg38_gencode_v27.txt --out-dir $D/haplodose --dataset $DS --code $CODE \
  --column log_ratio_gene --ref-holdout 0.3 --seed 20261103 --hd-args "$HD" 2> $D/run_haplodose_time.txt
$PY -s -m pip freeze > $D/pip_freeze.txt
echo -e "end\t$(date -Is)" >> $D/job.tsv
echo HAPLODOSE_DONE
