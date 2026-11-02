#!/usr/bin/env bash
#SBATCH --job-name=haplodose_comparators_smoke2
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/comparators_smoke2_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/comparators_smoke2_%j.err
# Second smoke: Haplodose-layout outputs for every runner, R inferCNV added, and the two scorers
# (cell_lines/evaluate.py, ucan_spatial/dna_truth.py) run on comparator outputs.
# Every tool run is its own srun step, so sacct reports per-step Elapsed and MaxRSS.
set -uo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME
# R/JAGS modules set PYTHONPATH to a system Python; keep only what Rscript needs, then purge.
# R's library path from a clean login shell: the batch env can carry other GCC libs first.
L=$(env -i HOME=$HOME bash -lc 'module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1; echo "$(command -v Rscript)|$LD_LIBRARY_PATH"')
# R's GCC 13.2 libstdc++ first, GCC 11/12 entries dropped, JAGS (rjags) last: older libstdc++ must not shadow R's GCC 13 one.
export RSCRIPT=${L%%|*} R_LD_LIBRARY_PATH=/nfs/sw/easybuild/software/GCCcore/13.2.0/lib64:$(echo "${L#*|}" | tr ':' '\n' | grep -v 'GCCcore/1[12]\.' | paste -sd:):/nfs/sw/easybuild/software/JAGS/4.3.2-foss-2023a/lib
module purge; unset PYTHONPATH PYTHONHOME
D=benchmarks/comparators
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
PROTO=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
BENCH=/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark
OUT=results/comparators_smoke2_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}
mkdir "$OUT"
printf 'step\texit\n' > "$OUT/status.tsv"
step() { local name=$1; shift; echo "== $name $(date)"
  srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --job-name="$name" "$@"
  printf '%s\t%s\n' "$name" "$?" >> "$OUT/status.tsv"; }
last() { ls -d "$1"/$2_* 2>/dev/null | tail -1; }

step simulate "$HD" -B -m haplodose.cli simulate --out "$OUT/sim"
SIM=(--counts "$OUT/sim/counts.h5" --cells "$OUT/sim/cells.tsv" --positions "$OUT/sim/positions.tsv")
for t in infercnvpy prototype infercnv_r; do
  py=$CNV; [[ $t == prototype ]] && py=$PROTO
  step sim_$t "$py" -B $D/run_$t.py "${SIM[@]}" --out-parent "$OUT/sim_runs"
done

SMC=(--counts "$BENCH/SMC_smoke/input/counts.h5" --cells "$BENCH/SMC_smoke/input/cells.tsv"
     --positions "$BENCH/SMC_smoke/input/positions.tsv")
step smc_infercnv_r "$CNV" -B $D/run_infercnv_r.py "${SMC[@]}" --exclude-genes /gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt --out-parent "$OUT/SMC_smoke_runs"

CL=results/celllines_20261101_0853/input
CLA=(--counts "$CL/counts.h5" --cells "$CL/cells.tsv" --positions "$CL/positions.tsv")
for t in infercnvpy prototype infercnv_r; do
  py=$CNV; [[ $t == prototype ]] && py=$PROTO
  step celllines_$t "$py" -B $D/run_$t.py "${CLA[@]}" --out-parent "$OUT/celllines_runs"
  r=$(last "$OUT/celllines_runs" "$t"); [[ -n $r ]] && step celllines_${t}_eval "$HD" -s -B benchmarks/cell_lines/evaluate.py "$r" "$CL" "$r/eval"
done

S=U0247; UI=$BENCH/$S/input
UA=(--counts "$UI/counts.h5" --cells "$UI/cells.tsv" --positions "$UI/positions.tsv")
for t in infercnvpy prototype; do
  py=$CNV; ex=(--exclude-genes /gpfs/commons/home/jameslee/CRC/src/44_haplodose/masks.txt)
  [[ $t == prototype ]] && { py=$PROTO; ex=(); }      # the prototype masks IGF2 itself
  step ${S}_$t "$py" -B $D/run_$t.py "${UA[@]}" "${ex[@]}" --out-parent "$OUT/${S}_runs"
  r=$(last "$OUT/${S}_runs" "$t"); [[ -n $r ]] && step ${S}_${t}_dna "$HD" -s -B benchmarks/ucan_spatial/dna_truth.py "$r" "$S" "$r/dna"
done
cat "$OUT/status.tsv"
sacct -j "$SLURM_JOB_ID" -o JobID%20,JobName%28,AllocCPUS,Elapsed,TotalCPU,MaxRSS,State -P > "$OUT/sacct.tsv" || true
