#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/comparators/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/comparators/%x_%j.err
# One comparator run, optionally followed by one scorer.
# Usage: sbatch -J cmp_<tool>_<dataset> [--mem=..] run_one.sbatch TOOL OUT_PARENT SCORER(none|cell_lines|ucan:<S>|moorman) -- RUNNER_ARGS...
# The tool runs as srun step 0, so `sacct -j <id>` gives its Elapsed, MaxRSS and AllocCPUS.
set -uo pipefail
TOOL=$1 PARENT=$2 SCORER=$3; shift 3; [[ $1 == -- ]] && shift
D=benchmarks/comparators
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
PROTO=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
if [[ $TOOL == infercnv_r ]]; then
  # R/JAGS modules set PYTHONPATH to a system Python; keep only what Rscript needs, then purge.
  # R's library path from a clean login shell: the batch env can carry other GCC libs first.
  L=$(env -i HOME=$HOME bash -lc 'module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1; echo "$(command -v Rscript)|$LD_LIBRARY_PATH"')
  # R's GCC 13.2 libstdc++ first, GCC 11/12 entries dropped, JAGS (rjags) last: older libstdc++ must not shadow R's GCC 13 one.
  export RSCRIPT=${L%%|*} R_LD_LIBRARY_PATH=/nfs/sw/easybuild/software/GCCcore/13.2.0/lib64:$(echo "${L#*|}" | tr ':' '\n' | grep -v 'GCCcore/1[12]\.' | paste -sd:):/nfs/sw/easybuild/software/JAGS/4.3.2-foss-2023a/lib
  module purge
fi
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=$CNV; [[ $TOOL == prototype ]] && PY=$PROTO
mkdir -p "$PARENT"
before=$(ls -d "$PARENT/${TOOL}_"* 2>/dev/null)
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --job-name="$TOOL" "$PY" -B $D/run_$TOOL.py "$@" --out-parent "$PARENT"
rc=$?
run=$(comm -13 <(echo "$before") <(ls -d "$PARENT/${TOOL}_"* 2>/dev/null) | tail -1)
echo "tool exit $rc run $run"
[[ $rc -eq 0 && -n $run ]] || exit $rc
echo "$SLURM_JOB_ID" > "$run/slurm_job_id"
case $SCORER in
  none) ;;
  cell_lines) $HD -s -B benchmarks/cell_lines/evaluate.py "$run" results/celllines_20261101_0853/input "$run/eval" ;;
  ucan:*) $HD -s -B benchmarks/ucan_spatial/dna_truth.py "$run" "${SCORER#ucan:}" "$run/dna" ;;
  moorman) $HD -s -B benchmarks/moorman/dna_truth.py "$run" "$run/dna" ;;
esac
