#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --job-name=masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 4 SPECTRUM: CellRanger count for one scRNA sample (all its SRA runs, each linked as one lane), to make the
# count matrix (input contract) and the BAM (Numbat allele pileup). Same reference, version and flags as the Schmid
# runs (benchmarks/schmid/cellranger/run_cellranger.sbatch, results/schmid_cellranger_20261102_2129): CellRanger
# 8.0.1, refdata-gex-GRCh38-2020-A, --create-bam=true, chemistry auto-detected, introns at the default.
# QC only: the job records cell calls and mapping metrics; no expression is inspected.
# Requires every run's FASTQ .ok marker (md5 checked by benchmarks/heldout2/fetch_spectrum.sbatch). Job name is
# masld-test (nslab convention); the id-to-sample map is in OUT/jobs.tsv and dev/job_renames_20261102.tsv.
# Usage: sbatch --output=OUT/logs/<sample>_%j.out run_cellranger.sbatch OUT SAMPLE_TITLE SAMPLES_TSV
set -euo pipefail
OUT=$1 TITLE=$2 TABLE=$3
H=/gpfs/commons/home/jameslee/Haplodose
FQ=$H/data/external/spectrum_ov_scrna/fastq
REF=/gpfs/commons/groups/sanjana_lab/jmorris/211215_BeeSTINGseq/downloads/refdata-gex-GRCh38-2020-A
row=$(awk -F'\t' -v t="$TITLE" '$1==t' "$TABLE")
[[ -n $row ]] || { echo "sample not in $TABLE: $TITLE"; exit 2; }
GROUP=$(cut -f3 <<< "$row"); IFS=, read -ra RUNS <<< "$(cut -f4 <<< "$row")"
SAMPLE=${RUNS[0]} SN=S1
test ! -e "$OUT/$TITLE" || { echo "$OUT/$TITLE exists; refusing to overwrite"; exit 1; }
S=$OUT/staged_fastq/$TITLE
mkdir -p "$S" "$OUT/logs"
lane=0
for r in "${RUNS[@]}"; do
  lane=$((lane + 1)); L=$(printf 'L%03d' $lane)
  for pair in 1:R1 2:R2; do
    i=${pair%%:*} rd=${pair##*:}
    f=$FQ/$GROUP/$r/${r}_$i.fastq.gz
    test -s "$f" && test -e "$f.ok" || { echo "MISSING or incomplete FASTQ $f"; exit 2; }
    ln -s "$f" "$S/${SAMPLE}_${SN}_${L}_${rd}_001.fastq.gz"
  done
done
module load CellRanger/8.0.1
cellranger --version > "$OUT/logs/${TITLE}_cellranger_version.txt"
md5sum "$REF/reference.json" > "$OUT/logs/${TITLE}_reference_md5.txt"
/usr/bin/time -v cellranger count --id="$TITLE" --output-dir="$OUT/$TITLE" \
  --transcriptome="$REF" --fastqs="$S" --sample="$SAMPLE" \
  --create-bam=true --localcores=16 --localmem=120 --disable-ui \
  2> "$OUT/logs/${TITLE}_time.txt"
O=$OUT/$TITLE/outs
test -s "$O/possorted_genome_bam.bam" && test -s "$O/possorted_genome_bam.bam.bai" && test -s "$O/raw_feature_bc_matrix.h5"
md5sum "$O/possorted_genome_bam.bam" "$O/raw_feature_bc_matrix.h5" "$O/filtered_feature_bc_matrix.h5" > "$OUT/$TITLE/outs_md5.txt"
touch "$OUT/$TITLE/done.ok"
