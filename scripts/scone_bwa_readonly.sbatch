#!/usr/bin/env bash
set -euo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Use an allocated compute node.' >&2; exit 2; }
TASK=$(realpath "$1")
source /etc/profile.d/modules.sh
module purge
module load BWA/0.7.18-GCCcore-12.3.0
export HAPLODOSE_BWA_BIN
HAPLODOSE_BWA_BIN=$(command -v bwa)
export HAPLODOSE_BWA_LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}
export HAPLODOSE_BWA_MODULE_RECORD=$TASK/bwa_modules_$SLURM_JOB_ID.txt
module -t list > "$HAPLODOSE_BWA_MODULE_RECORD" 2>&1
exec /usr/bin/bash "$TASK/analysis_readonly.sbatch" "$@"
