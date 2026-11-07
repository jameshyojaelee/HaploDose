#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=b6k
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/scaling/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/scaling/%x_%j.err
# GPU arm of the scaling benchmark: HaploDose only, same inputs and outputs layout as run_one.sbatch, written to
# ROOT/runs/haplodose_gpu_<N>_rep<k>. Hard --constraint=b6k (not --prefer): GPU timings are compared across sizes,
# so every run must use the same GPU model (RTX PRO 6000 Blackwell; smoke 22612305 passed there). 16 CPUs to match
# the CPU arm (Leiden, PCA and I/O stay on the CPU). JAX_PLATFORMS=cuda: the frozen code predates --device and
# binds to whatever JAX_PLATFORMS names. Device memory is sampled every second by nvidia-smi into gpu_mem.tsv
# (MiB, whole GPU; preallocation is off so it tracks use). Host memory is set per job on the sbatch line.
# Usage: sbatch -J scl_haplodose_gpu_<N>_<rep> --mem=<M> run_gpu.sbatch ROOT N REP
set -euo pipefail
ROOT=$1 N=$2 REP=$3
I=$ROOT/input O=$ROOT/runs/haplodose_gpu_${N}_rep${REP}
test ! -e "$O"; mkdir -p "$O"
HD=${HAPLODOSE_PYTHON:-/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python}
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=16
export PYTHONPATH=${HAPLODOSE_CODE:-$ROOT/code} JAX_PLATFORMS=cuda XLA_PYTHON_CLIENT_PREALLOCATE=false
if [[ ${HAPLODOSE_WHEEL_ONLY:-0} == 1 ]]; then unset PYTHONPATH; fi
export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
echo -e "host\t$(hostname)\nslurm_job\t$SLURM_JOB_ID\ncpus\t$SLURM_CPUS_PER_TASK\ngpu\t$GPU\nstart\t$(date -Is)" > "$O/job.tsv"
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv,noheader,nounits -l 1 > "$O/gpu_mem.tsv" &
SMI=$!
record_exit() {
  status=$?
  kill "$SMI" 2>/dev/null || true
  printf 'exit\t%s\nend\t%s\n' "$status" "$(date -Is)" >> "$O/job.tsv"
  trap - EXIT
  exit "$status"
}
trap record_exit EXIT
HELP=$($HD -s -m haplodose.cli run --help)
EXTRA=()
if [[ $HELP == *--device* ]]; then EXTRA+=(--device gpu); fi
if [[ $HELP == *--counts-kind* ]]; then
  EXTRA+=(--counts-kind raw --genome-build hg38 --coordinate-system one-based-closed --gene-namespace custom)
fi
/usr/bin/time -v $HD -s -m haplodose.cli run --counts "$I/counts_$N.h5" --cells "$I/cells_$N.tsv" \
  --positions "$I/positions.tsv" --out "$O/out" "${EXTRA[@]}" 2> "$O/time.txt"
$HD -s -c 'import json,sys; q=json.load(open(sys.argv[1])); assert q["backend"] == "gpu", q["backend"]' "$O/out/qc.json"
