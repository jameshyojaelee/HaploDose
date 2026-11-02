#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/andor/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/andor/%x_%j.err
# Andor 2020 10x Single Cell CNV (PRJNA498809): fetch one run and count reads per cell in 500 kb bins.
# The SRA run is the submitters' Cell Ranger DNA BAM loaded as an alignment database (bam-load), so
# alignments (GRCh38, chr names) are reused as is; no realignment. Cells are read groups: each @RG has
# the cell barcode as SM (checked on SNU-638 and in the original SNU-601 BAM).
# Read filter: primary, read 1 of a proper pair (flag 0x42 set), not unmapped/secondary/QC-fail/
# duplicate/supplementary (0xF04 clear), MAPQ >= 30. One count per fragment.
# Usage: sbatch -J andor_fetch_<line> fetch_bins.sbatch LINE RUN
set -euo pipefail
LINE=$1 RUN=$2
D=data/external/andor2020_scdna
mkdir -p $D/sra $D/bins
module load SRA-Toolkit/3.2.1 samtools/1.23.1
# Resumes; also fetches reference contigs an earlier run missed. Up to 3 tries; a contig that still fails
# (MKN-45: unplaced KI/GL contigs, "cannot get remote location") is left out, since only chr1-22 are counted.
for i in 1 2 3; do prefetch --max-size 200G -O $D/sra $RUN && break || { echo "prefetch try $i failed"; sleep 60; }; done
vdb-validate $D/sra/$RUN/$RUN.sra
OUT=$D/bins/${LINE}.tsv.gz
[[ ! -e $OUT ]] || { echo "exists: $OUT"; exit 0; }
# --header-only printed nothing on a remote run; a tiny region query prints the full header.
sam-dump --aligned-region chr1:1000000-1001000 $D/sra/$RUN/$RUN.sra | grep '^@RG' > $D/bins/${LINE}.rg.txt
sam-dump --primary $D/sra/$RUN/$RUN.sra \
  | samtools view -f 0x42 -F 0xF04 -q 30 - \
  | awk -F'\t' -v OFS='\t' '{rg=""; for (i=12; i<=NF; i++) if (substr($i,1,5)=="RG:Z:") {rg=substr($i,6); break}
        if (rg=="") next; c[rg OFS $3 OFS int($4/500000)]++; n++}
        END {for (k in c) print k, c[k]; print "#fragments", n > "/dev/stderr"}' \
  | gzip -c > $OUT.tmp
mv $OUT.tmp $OUT
echo "done $LINE $RUN $(zcat $OUT | wc -l) cell-bin rows"
