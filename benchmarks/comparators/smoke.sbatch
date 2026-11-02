#!/usr/bin/env bash
#SBATCH --job-name=haplodose_comparators_smoke
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/comparators_smoke_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/comparators_smoke_%j.err
# Smoke test of the comparator runners: simulator output (infercnvpy, prototype),
# SMC_smoke (infercnvpy, prototype), YCC061 with alleles (Numbat). Each step runs even
# if an earlier one fails; status.tsv records exit codes.
set -uo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME
D=benchmarks/comparators
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
PROTO=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
BENCH=/gpfs/commons/home/jameslee/CRC/results/cnv_caller_benchmark
OUT=results/comparators_smoke_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}
mkdir "$OUT"
printf 'step\texit\n' > "$OUT/status.tsv"
step() { local name=$1; shift; echo "== $name $(date)"; "$@"; printf '%s\t%s\n' "$name" "$?" >> "$OUT/status.tsv"; }

step simulate "$HD" -B -m haplodose.cli simulate --out "$OUT/sim"
SIM=(--counts "$OUT/sim/counts.h5" --cells "$OUT/sim/cells.tsv" --positions "$OUT/sim/positions.tsv")
step sim_infercnvpy "$CNV" -B $D/run_infercnvpy.py "${SIM[@]}" --out-parent "$OUT/sim_runs"
step sim_prototype "$PROTO" -B $D/run_prototype.py "${SIM[@]}" --out-parent "$OUT/sim_runs"

SMC=(--counts "$BENCH/SMC_smoke/input/counts.h5" --cells "$BENCH/SMC_smoke/input/cells.tsv"
     --positions "$BENCH/SMC_smoke/input/positions.tsv")
printf 'IGF2\n' > "$OUT/exclude_igf2.txt"
step smc_infercnvpy "$CNV" -B $D/run_infercnvpy.py "${SMC[@]}" --exclude-genes "$OUT/exclude_igf2.txt" --out-parent "$OUT/SMC_smoke_runs"
step smc_prototype "$PROTO" -B $D/run_prototype.py "${SMC[@]}" --out-parent "$OUT/SMC_smoke_runs"

# Imprinted-gene exclusion copied from T2_run.R.
printf '%s\n' IGF2 H19 KCNQ1OT1 CDKN1C KCNQ1 PHLDA2 SLC22A18 MEG3 DLK1 MEST PEG3 PEG10 SNRPN SNURF NNAT PLAGL1 GNAS > "$OUT/exclude_imprinted.txt"
step ycc061_numbat "$CNV" -B $D/run_numbat.py \
    --counts /gpfs/commons/home/jameslee/CRC/data/external/gse273567/YCC061/alignment/YCC061/outs/filtered_feature_bc_matrix.h5 \
    --cells "$BENCH/YCC061/input/cells.tsv" --positions "$BENCH/YCC061/input/positions.tsv" \
    --alleles /gpfs/commons/home/jameslee/CRC/results/numbat_11p15/YCC061/alleles/YCC061_allele_counts.tsv.gz \
    --specimen YCC061 --exclude-genes "$OUT/exclude_imprinted.txt" --ncores 8 --out-parent "$OUT/YCC061_runs"
cat "$OUT/status.tsv"
