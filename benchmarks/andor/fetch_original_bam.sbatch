#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/andor/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/andor/%x_%j.err
# Per-cell 500 kb bin counts from the submitters' original Cell Ranger DNA BAM (public only for SNU-601,
# bam_sources.tsv), streamed from S3 without storing it. Here every read has its cell read group (RG:Z:<id>,
# header SM = barcode) and CB tag. The SRA alignment load keeps the cell read group only for a few
# thousand reads per cell; the rest carry the 8-mer sample index as spot group (checked on SNU-638).
# Same read filter as fetch_bins.sbatch: flag 0x42 set, 0xF04 clear, MAPQ >= 30, one count per fragment.
# Usage: sbatch -J andor_bam_<line> fetch_original_bam.sbatch LINE URL MD5
set -euo pipefail
LINE=$1 URL=$2 MD5=$3
D=data/external/andor2020_scdna/bins_original_bam
mkdir -p $D
OUT=$D/${LINE}.tsv.gz
[[ ! -e $OUT ]] || { echo "exists: $OUT"; exit 1; }
module load samtools/1.23.1
curl -sS --retry 5 "$URL" \
  | tee >(md5sum | cut -d' ' -f1 > $D/${LINE}.md5) \
  | samtools view -h -f 0x42 -F 0xF04 -q 30 - \
  | awk -F'\t' -v OFS='\t' -v H=$D/${LINE}.rg.txt '/^@RG/ {print > H; next} /^@/ {next}
        {rg=""; for (i=12; i<=NF; i++) if (substr($i,1,5)=="RG:Z:") {rg=substr($i,6); break}
         if (rg=="") next; c[rg OFS $3 OFS int($4/500000)]++; n++}
        END {for (k in c) print k, c[k]; print "#fragments", n > "/dev/stderr"}' \
  | gzip -c > $OUT.tmp
sleep 5
got=$(cat $D/${LINE}.md5)
[[ $got == "$MD5" ]] || { echo "md5 mismatch: $got vs $MD5"; exit 1; }
mv $OUT.tmp $OUT
echo "done $LINE md5 ok $(zcat $OUT | wc -l) cell-bin rows"
