#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --constraint=v6
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/scaling/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/scaling/%x_%j.err
# Compute scaling benchmark, one tool x size x repeat. Hard --constraint=v6 (not --prefer): timings are compared
# across tools, so every run must be on the same CPU generation; R inferCNV also needs it (SIGILL on v4).
# Memory is set per job on the sbatch line. 72 h cap: a run that hits it is reported as not finished.
# Usage: sbatch -J scl_<tool>_<N>_<rep> --mem=<M> run_one.sbatch ROOT TOOL N REP
set -euo pipefail
ROOT=$1 TOOL=$2 N=$3 REP=$4
I=$ROOT/input O=$ROOT/runs/${TOOL}_${N}_rep${REP}
test ! -e "$O"; mkdir -p "$O"
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=16
echo -e "host\t$(hostname)\nslurm_job\t$SLURM_JOB_ID\ncpus\t$SLURM_CPUS_PER_TASK\nstart\t$(date -Is)" > "$O/job.tsv"
ARGS=(--counts "$I/counts_$N.h5" --cells "$I/cells_$N.tsv" --positions "$I/positions.tsv")
set +e
case $TOOL in
haplodose)
  export PYTHONPATH=$ROOT/code JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
  export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
  /usr/bin/time -v $HD -s -m haplodose.cli run "${ARGS[@]}" --out "$O/out" 2> "$O/time.txt" ;;
infercnvpy)
  /usr/bin/time -v $CNV -s benchmarks/comparators/run_infercnvpy.py "${ARGS[@]}" --out-parent "$O" 2> "$O/time.txt" ;;
infercnv_r)
  L=$(env -i HOME=$HOME bash -lc 'module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1; echo "$(command -v Rscript)|$LD_LIBRARY_PATH"')
  export RSCRIPT=${L%%|*} R_LD_LIBRARY_PATH=/nfs/sw/easybuild/software/GCCcore/13.2.0/lib64:$(echo "${L#*|}" | tr ':' '\n' | grep -v 'GCCcore/1[12]\.' | paste -sd:):/nfs/sw/easybuild/software/JAGS/4.3.2-foss-2023a/lib
  unset PYTHONPATH PYTHONHOME
  /usr/bin/time -v $CNV -s benchmarks/comparators/run_infercnv_r.py "${ARGS[@]}" --out-parent "$O" 2> "$O/time.txt" ;;
copykat|scevan)
  unset PYTHONPATH PYTHONHOME R_LIBS R_LIBS_USER
  /usr/bin/time -v $CNV -s benchmarks/comparators/run_r_field.py --tool $TOOL "${ARGS[@]}" --out-parent "$O" 2> "$O/time.txt" ;;
esac
echo -e "exit\t$?\nend\t$(date -Is)" >> "$O/job.tsv"
