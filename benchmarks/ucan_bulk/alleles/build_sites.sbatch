#!/usr/bin/env bash
#SBATCH --job-name=ucan_allele_sites
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_allele_sites_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/ucan_allele_sites_%j.err
# Common exonic SNPs for RNA allele counting: 1000 Genomes phase 3 on GRCh38 (no chr prefix, same naming as the
# U-CAN RNA CRAMs), biallelic SNVs with EUR_AF 0.05-0.95, inside exons of Ensembl 112 protein-coding and lncRNA
# genes on chr1-22. Usage: sbatch build_sites.sbatch OUTDIR
set -euo pipefail
OUT=$1
mkdir -p "$OUT"
module load BCFtools/1.19-GCC-13.2.0
G=/gpfs/commons/datasets/1000genomes/hg38
GTF=/gpfs/commons/home/jameslee/CRC/data/external/ensembl112/Homo_sapiens.GRCh38.112.gtf.gz
zcat $GTF | awk -F'\t' '$3=="exon" && $1 ~ /^([1-9]|1[0-9]|2[0-2])$/ && ($9 ~ /gene_biotype "protein_coding"/ || $9 ~ /gene_biotype "lncRNA"/) {print $1"\t"$4-1"\t"$5}' \
  | sort -k1,1V -k2,2n | awk 'BEGIN{OFS="\t"} $1==c && $2<=e {if($3>e)e=$3; next} {if(c!="")print c,s,e; c=$1;s=$2;e=$3} END{print c,s,e}' > "$OUT/exons.bed"
one() {
  bcftools view -v snps -m2 -M2 -i 'EUR_AF>=0.05 && EUR_AF<=0.95' -R "$OUT/exons.bed" \
    "$G/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz" -Ou \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%EUR_AF\n' > "$OUT/sites.chr$1.tsv"
}
export -f one; export OUT G
seq 1 22 | xargs -P 8 -I{} bash -c 'one {}'
for c in $(seq 1 22); do cat "$OUT/sites.chr$c.tsv"; done | sort -k1,1V -k2,2n -u > "$OUT/sites.tsv"
rm "$OUT"/sites.chr*.tsv
cut -f1,2 "$OUT/sites.tsv" > "$OUT/sites.regions.tsv"
awk 'BEGIN{OFS="\t"}{print $1,$2,$4","$5}' "$OUT/sites.tsv" | bgzip > "$OUT/sites.alleles.tsv.gz"
tabix -s1 -b2 -e2 "$OUT/sites.alleles.tsv.gz"
wc -l "$OUT/exons.bed" "$OUT/sites.tsv"
bcftools --version | head -1 > "$OUT/versions.txt"
