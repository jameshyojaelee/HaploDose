#!/usr/bin/env bash
set -euo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Use an allocated compute node.' >&2; exit 2; }
module purge
module load R/4.3.3-gfbf-2023b
module -t list 2> "$1/r_modules_$SLURM_JOB_ID.txt"
exec /usr/bin/bash /gpfs/commons/home/jameslee/Haplodose/dev/scone_count_chain_20261109/analysis_readonly.sbatch "$@"
