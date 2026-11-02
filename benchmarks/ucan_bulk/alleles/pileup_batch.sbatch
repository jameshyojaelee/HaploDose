#!/usr/bin/env bash
#SBATCH --job-name=ucan_pileup
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/%x_%j.out
# Run pileup.sh (MODE T, streaming) for every sample in SAMPLE_LIST, 7 at a time (2 threads each).
# Usage: sbatch pileup_batch.sbatch SAMPLE_LIST SITES_DIR OUTDIR
set -euo pipefail
LIST=$1 SITES=$2 OUT=$3
source /etc/profile.d/modules.sh; module load BCFtools/1.19-GCC-13.2.0
P=/gpfs/commons/home/jameslee/Haplodose/benchmarks/ucan_bulk/alleles/pileup.sh
xargs -a "$LIST" -P 7 -I{} bash -c "MODE=T bash $P {} $SITES $OUT || echo FAILED {}"
echo "done $(ls $OUT/*.alleles.tsv.gz | wc -l) files in $OUT"
