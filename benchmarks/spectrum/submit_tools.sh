#!/usr/bin/env bash
# Claim 4 SPECTRUM: submit the allele pileups, then every tool per patient (after its pileup for numbat and
# haplodose_alleles), then the paired test after all tool jobs. One job per patient and tool, no arrays.
# Run only after the lead has written the reference-cell rule into docs/STATUS.md and prepare_inputs.sbatch is done.
# Job ids go to ROOT/jobs.tsv and dev/job_renames_20261102.tsv (jobs run as masld-test on nslab).
# Usage: submit_tools.sh ROOT CR_ROOT
set -euo pipefail
ROOT=$1 CR=$2
H=/gpfs/commons/home/jameslee/Haplodose S=$H/benchmarks/spectrum MAP=$H/dev/job_renames_20261102.tsv
mkdir -p "$ROOT/logs"
[[ -f $ROOT/jobs.tsv ]] || echo -e "jobid\ttask\tsubmitted" > "$ROOT/jobs.tsv"
rec() { echo -e "$1\t$2\t$(date -Is)" >> "$ROOT/jobs.tsv"; echo -e "$1\tsp_$2\tmasld-test\tnslab" >> "$MAP"; echo "$1 $2"; }
declare -A RES=([haplodose]="--cpus-per-task=8 --mem=64G" [haplodose_alleles]="--cpus-per-task=8 --mem=64G"
                [infercnvpy]="--cpus-per-task=8 --mem=64G" [copykat]="--cpus-per-task=16 --mem=128G"
                [scevan]="--cpus-per-task=16 --mem=128G" [numbat]="--cpus-per-task=16 --mem=96G")
ALL=()
for d in "$ROOT"/input/SPECTRUM-OV-*/; do
  P=$(basename "$d")
  pu=$(sbatch --parsable --output="$ROOT/logs/${P}_pileup_%j.out" $S/pileup_phase.sbatch "$ROOT/input" "$CR" "$ROOT/alleles" "$P")
  rec "$pu" "${P}_pileup"
  for T in haplodose infercnvpy copykat scevan numbat haplodose_alleles; do
    DEP=(); [[ $T == numbat || $T == haplodose_alleles ]] && DEP=(--dependency=afterok:$pu)
    id=$(sbatch --parsable ${RES[$T]} "${DEP[@]}" --output="$ROOT/logs/${P}_${T}_%j.out" $S/run_tool.sbatch "$ROOT" "$P" "$T")
    rec "$id" "${P}_$T"; ALL+=("$id")
  done
done
dep=$(IFS=:; echo "${ALL[*]}")
id=$(sbatch --parsable --account=nslab --qos=nslab -J masld-test --partition=cpu --prefer=v6 --cpus-per-task=1 --mem=8G --time=2:00:00 \
  --chdir=$H --dependency=afterany:$dep --output="$ROOT/logs/paired_test_%j.out" \
  --wrap "/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s benchmarks/spectrum/paired_test.py $H/dev/spectrum_run_plan_20261103_0750/patients_with_scrna.tsv $ROOT/scores $ROOT/primary")
rec "$id" paired_test
