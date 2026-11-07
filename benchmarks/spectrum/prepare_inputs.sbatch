#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --job-name=masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 4 SPECTRUM: per-patient input contract (prepare_inputs.py) once every CellRanger sample is done.
# Usage: sbatch --output=ROOT/logs/prepare_%j.out prepare_inputs.sbatch ROOT CR_ROOT
set -euo pipefail
ROOT=$1 CR=$2
PLAN=dev/spectrum_run_plan_20261103_0750
GEO=data/external/spectrum_ov_scrna/geo/GSE180661_GEO_cells.tsv.gz
sha256sum --quiet -c $PLAN/code_sha256.txt || { echo "eval code hash mismatch" >&2; exit 2; }
[[ $(md5sum < $GEO) == "3725f8259e0bfdff61632c2735ad7286  -" ]] || { echo "GEO cell table md5 differs" >&2; exit 2; }
for s in $(awk -F'\t' 'NR>1{print $1}' $PLAN/samples.tsv); do test -e "$CR/$s/done.ok" || { echo "not done: $s" >&2; exit 1; }; done
/usr/bin/time -v /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s benchmarks/spectrum/prepare_inputs.py \
  $PLAN/samples.tsv $GEO "$CR" "$ROOT/input"
