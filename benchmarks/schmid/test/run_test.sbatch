#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# End-to-end test of the Schmid-pipeline HaploDose wrapper on Kinker 2020 development data (not Schmid data):
# build their input layout (pool18, rotation r0, query line C32_SKIN), write count_matrix.txt with R write.table,
# run benchmarks/schmid/run_haplodose.py on the frozen code (CPU), then score with their lib.R + lib_haplodose.R.
# Usage: sbatch --output=OUT/job_%j.out benchmarks/schmid/test/run_test.sbatch OUT
set -euo pipefail
OUT=$1 DS=KINKERC32
H=/gpfs/commons/home/jameslee/Haplodose
F=$H/results/heldout_20261102_1424
CODE=$H/results/scaling_frozen1424_20261102_1541/code
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
T=$H/benchmarks/schmid/test
(cd $F && sha256sum --quiet -c code_sha256.txt)
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ncpu\t$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $OUT/job.tsv
$PY -I $T/build_kinker_input.py $OUT pool18 C32_SKIN $DS
/usr/bin/time -v $R $T/write_count_matrix.R $OUT/input_$DS 2> $OUT/write_count_matrix_time.txt
rm $OUT/input_$DS/counts.mtx
/usr/bin/time -v $PY -I $H/benchmarks/schmid/run_haplodose.py \
  --matrix $OUT/input_$DS/count_matrix.txt --annot $OUT/input_$DS/sample_annotation.txt \
  --ref-groups $OUT/input_$DS/ref_groups.txt --gene-pos $H/external/annotations/hg38_gencode_v27.txt \
  --out-dir $OUT/output/$DS/haplodose --dataset $DS --code $CODE 2> $OUT/run_haplodose_time.txt
$R $T/score_test.R $OUT $DS
$PY -s -m pip freeze > $OUT/pip_freeze.txt
echo TEST_DONE
