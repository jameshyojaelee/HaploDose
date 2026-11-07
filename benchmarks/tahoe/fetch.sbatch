#!/usr/bin/env bash
#SBATCH --job-name=tahoe_fetch
#SBATCH --account=innovation
#SBATCH --qos=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/tahoe/%x_%j.out
# Tahoe-100M (tahoebio/Tahoe-100M, CC0) metadata, plate1 DMSO_TF cell selection and the expression shards that hold
# those cells (benchmarks/tahoe/fetch_tahoe.py; claim 8, docs/STATUS.md). Resumable: rerunning skips finished steps.
set -euo pipefail
D=data/external/tahoe100m_20261106
mkdir -p $D
PY=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
$PY -I -m pip freeze > $D/pip_freeze_fetch.txt
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\nstart\t$(date -Is)" > $D/fetch_job.tsv
$PY -I benchmarks/tahoe/fetch_tahoe.py $D
sha256sum benchmarks/tahoe/fetch_tahoe.py >> $D/fetch_job.tsv
du -sb $D/metadata $D/shards $D/plate1_dmso_tf.parquet >> $D/fetch_job.tsv
echo -e "end\t$(date -Is)" >> $D/fetch_job.tsv
