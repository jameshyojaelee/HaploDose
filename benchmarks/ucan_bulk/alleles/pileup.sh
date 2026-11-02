#!/usr/bin/env bash
# Allele counts at common exonic SNPs in one U-CAN RNA alignment (BAM or CRAM, read only).
# Output TSV.gz: chrom, pos, ref, alt (1000G), ref_count, alt_count, depth (all bases at the site).
# bcftools mpileup: MAPQ >= 20, base quality >= 20, no BAQ (RNA), indels skipped, orphan pairs kept, depth cap 2000
# (counted once per read pair by mpileup's overlap handling).
# Usage: pileup.sh SAMPLE SITES_DIR OUTDIR
set -euo pipefail
command -v bcftools >/dev/null || { source /etc/profile.d/modules.sh; module load BCFtools/1.19-GCC-13.2.0; }
S=$1 SITES=$2 OUT=$3
CRC=/gpfs/commons/home/jameslee/CRC
# Indexed BAM where one exists, otherwise the CRAM with its index from interim/crai (htslib ##idx## syntax).
if test -f "$CRC/data/raw/bam/$S.bam.bai"; then ALN=$CRC/data/raw/bam/$S.bam
else ALN="$CRC/data/raw/cram/$S.cram##idx##$CRC/data/interim/crai/$S.cram.crai"; fi
REF=/gpfs/commons/home/jameslee/CRC/data/external/ensembl112/Homo_sapiens.GRCh38.dna.primary_assembly.fa
test -f "$OUT/$S.alleles.tsv.gz" && exit 0
# MODE R: index jumps to each site; MODE T (default): stream the file, keep only target sites. Depth capped at 2000.
if test "${MODE:-T}" = R; then SEL=(-R "$SITES/sites.regions.tsv"); else SEL=(-T "$SITES/sites.regions.tsv"); fi
bcftools mpileup -f "$REF" "${SEL[@]}" -a AD,DP -q 20 -Q 20 -B -I -A -d 2000 --threads 2 \
    -Ou "$ALN" 2> "$OUT/$S.mpileup.log" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' \
  | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {alt[$1":"$2]=$5; next}
       ($1":"$2) in alt { a=alt[$1":"$2]; n=split($4,al,","); split($5,ad,","); ac=0; tot=0;
         for(i=1;i<=n+1;i++) tot+=ad[i]; for(i=1;i<=n;i++) if(al[i]==a) ac=ad[i+1];
         if (tot>0) print $1,$2,$3,a,ad[1],ac,tot }' "$SITES/sites.tsv" - \
  | gzip > "$OUT/$S.alleles.tsv.gz.tmp"
mv "$OUT/$S.alleles.tsv.gz.tmp" "$OUT/$S.alleles.tsv.gz"
