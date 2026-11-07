#!/usr/bin/env bash
# Submit run_cellranger.sbatch for every sample in SAMPLES_TSV whose patient's FASTQs (all samples, all runs) have
# .ok markers and that has not been submitted yet (OUT/jobs.tsv). Safe to rerun as the download progresses. One job per sample, no arrays.
# Records each job id in OUT/jobs.tsv and dev/job_renames_20261102.tsv (jobs run as masld-test on nslab).
# Usage: submit_cellranger.sh OUT SAMPLES_TSV
set -euo pipefail
OUT=$1 TABLE=$2
H=/gpfs/commons/home/jameslee/Haplodose
FQ=$H/data/external/spectrum_ov_scrna/fastq
MAP=$H/dev/job_renames_20261102.tsv
mkdir -p "$OUT/logs"
[[ -f $OUT/jobs.tsv ]] || echo -e "jobid\tsample_title\tsubmitted" > "$OUT/jobs.tsv"
# Patients with any FASTQ not yet verified (.ok missing) wait.
declare -A NOTREADY
while IFS=$'\t' read -r title patient group runs gb ncells; do
  [[ $title == sample_title ]] && continue
  for r in ${runs//,/ }; do
    for i in 1 2; do [[ -e $FQ/$group/$r/${r}_$i.fastq.gz.ok ]] || NOTREADY[$patient]=1; done
  done
done < "$TABLE"
n=0
while IFS=$'\t' read -r title patient group runs gb ncells; do
  [[ $title == sample_title ]] && continue
  grep -qP "\t$title\t" "$OUT/jobs.tsv" && continue
  [[ -z ${NOTREADY[$patient]:-} ]] || continue
  id=$(sbatch --parsable --output="$OUT/logs/${title}_%j.out" "$H/benchmarks/spectrum/run_cellranger.sbatch" "$OUT" "$title" "$TABLE")
  echo -e "$id\t$title\t$(date -Is)" >> "$OUT/jobs.tsv"
  echo -e "$id\tcr_$title\tmasld-test\tnslab" >> "$MAP"
  echo "submitted $id $title"; n=$((n + 1))
done < "$TABLE"
echo "patients waiting for FASTQs: ${#NOTREADY[@]}"
echo "$n submitted; $(($(wc -l < "$OUT/jobs.tsv") - 1)) of $(($(wc -l < "$TABLE") - 1)) samples submitted in total"
