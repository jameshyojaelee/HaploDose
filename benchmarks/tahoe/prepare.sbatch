#!/usr/bin/env bash
#SBATCH --job-name=tahoe_input
#SBATCH --account=innovation
#SBATCH --qos=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/tahoe/%x_%j.out
# Claim 8 Tahoe input (benchmarks/tahoe/prepare_tahoe.py) from the fetch_tahoe.py extraction into a new directory.
set -euo pipefail
SRC=data/external/tahoe100m_20261106
OUT=results/tahoe_input_$(date +%Y%m%d_%H%M)
test ! -e $OUT
test -f $SRC/verify.json
PY=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
$PY -I benchmarks/tahoe/prepare_tahoe.py $SRC $OUT
$PY -I -m pip freeze > $OUT/pip_freeze.txt
sha256sum benchmarks/tahoe/prepare_tahoe.py benchmarks/tahoe/fetch_tahoe.py benchmarks/kinker/prepare.py > $OUT/code_sha256.txt
cp $SRC/selection_rule.txt $SRC/pass_filter_values.tsv $SRC/verify.json $SRC/manifest.tsv $OUT/
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\nend\t$(date -Is)" > $OUT/job.tsv
