#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5: CellRanger count per Schmid dataset, to make the BAMs that CaSpER (BAFExtract) and Numbat (cellsnp-lite)
# need. PLANNED, NOT SUBMITTED: hold until the lead confirms (benchmarks/schmid/PLAN.md, "CellRanger plan").
# Mirrors external/benchmark_scrnaseq_cnv_callers/preprocessing_datasets/cellranger_*.sh (cellranger-7.0.0,
# refdata-gex-GRCh38-2020-A, --sample=<run>, COLO320 --chemistry=ARC-v1, localcores 6 localmem 64) with these
# recorded deviations: CellRanger 8.0.1 (STATUS claim 5 amendment f); --create-bam=true (required in 8.x; 7.0.0
# always wrote a BAM); 16 cores / 120 GB instead of 6 / 64 (speed only). --include-introns stays at the default
# (true in both 7.0.0 and 8.0.1). FASTQs are linked, never moved, under the bcl2fastq names the authors used.
# Touches only FASTQs; writes BAMs and matrices to OUT/<DS>/.
# Usage: sbatch -J cr_<DS> --output=OUT/logs/%x_%j.out run_cellranger.sbatch OUT DS
set -euo pipefail
OUT=$1 DS=$2
H=/gpfs/commons/home/jameslee/Haplodose
FQ=$H/data/external/schmid2025_raw/fastq
REF=/gpfs/commons/groups/sanjana_lab/jmorris/211215_BeeSTINGseq/downloads/refdata-gex-GRCh38-2020-A
CHEM=()
case $DS in
  COLO320)   RUNS=(SRR12900668); SAMPLE=SRR12900668; SN=S14; DIR=$FQ/COLO320HSR; CHEM=(--chemistry=ARC-v1) ;;
  BCC06)     RUNS=(SRR8315756);  SAMPLE=SRR8315756;  SN=S1;  DIR=$FQ/BCC06 ;;
  BCC06post) RUNS=(SRR8315759);  SAMPLE=SRR8315759;  SN=S1;  DIR=$FQ/BCC06post ;;
  MCF7)      RUNS=(SRR10018149); SAMPLE=SRR10018149; SN=S1;  DIR=$FQ/MCF7 ;;
  MM)        RUNS=(SRR10574405 SRR10574406 SRR10574407 SRR10574408 SRR10574409 SRR10574410 SRR10574411 SRR10574412)
             SAMPLE=RTBM199T; SN=S14; DIR=$FQ/MM199 ;;
  *) echo "unknown dataset $DS"; exit 2 ;;
esac
test ! -e "$OUT/$DS" || { echo "$OUT/$DS exists; refusing to overwrite"; exit 1; }
S=$OUT/staged_fastq/$DS
mkdir -p "$S" "$OUT/logs"
# Link each run as one lane (MM: 8 runs -> L001..L008 with I1, as in download_MM.sh). Require the fetch .ok marker.
lane=0
for r in "${RUNS[@]}"; do
  lane=$((lane + 1)); L=$(printf 'L%03d' $lane)
  reads=(1:R1 2:R2); [ "$DS" = MM ] && reads+=(3:I1)
  for pair in "${reads[@]}"; do
    i=${pair%%:*} rd=${pair##*:}
    if [ "$DS" = MM ]; then
      # MM199 fetch layout is flat, with one per-run marker $DIR/<run>.ok listing each validated file.
      f=$DIR/${r}_$i.fastq.gz
      test -s "$f" && grep -q "/${r}_$i.fastq.gz"$'\t' "$DIR/$r.ok" || { echo "MISSING or incomplete FASTQ $f"; exit 2; }
    else
      f=$DIR/$r/${r}_$i.fastq.gz
      test -s "$f" && test -e "$f.ok" || { echo "MISSING or incomplete FASTQ $f"; exit 2; }
    fi
    ln -s "$f" "$S/${SAMPLE}_${SN}_${L}_${rd}_001.fastq.gz"
  done
done
module load CellRanger/8.0.1
cellranger --version > "$OUT/logs/${DS}_cellranger_version.txt"
md5sum "$REF/reference.json" > "$OUT/logs/${DS}_reference_md5.txt"
/usr/bin/time -v cellranger count --id="${DS}_hg38_mapped" --output-dir="$OUT/$DS" \
  --transcriptome="$REF" --fastqs="$S" --sample="$SAMPLE" "${CHEM[@]}" \
  --create-bam=true --localcores=16 --localmem=120 --disable-ui \
  2> "$OUT/logs/${DS}_time.txt"
test -s "$OUT/$DS/outs/possorted_genome_bam.bam" && test -s "$OUT/$DS/outs/possorted_genome_bam.bam.bai"
md5sum "$OUT/$DS/outs/possorted_genome_bam.bam" "$OUT/$DS/outs/filtered_feature_bc_matrix.h5" > "$OUT/$DS/outs_md5.txt"
