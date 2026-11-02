#!/usr/bin/env bash
#SBATCH --job-name=numbat_hg19_panel
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Numbat hg19 phasing inputs, built because the authors' prebuilt hg19 panel URL is bot-walled:
# - panel/chr{1..22}.genotypes.bcf: 1000G phase 3 (20130502, hg19) biallelic SNPs, contigs renamed to chrN
#   (pileup_and_phase.R writes phasing targets with a chr prefix);
# - snps_AF5e2_chr1to22_hg19.vcf.gz: sites-only, biallelic SNPs with AF >= 0.05, contigs 1..22 (as the BAMs);
# - genetic_map_hg19_withX.txt.gz from Eagle 2.4.1 tables (its name tells pileup_and_phase.R genome = hg19).
set -euo pipefail
module load BCFtools/1.19-GCC-13.2.0
D=data/external/numbat_hg19_panel
mkdir -p $D/panel $D/raw $D/sites
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
[[ -e $D/genetic_map_hg19_withX.txt.gz ]] || curl -sS -L --retry 5 -o $D/genetic_map_hg19_withX.txt.gz \
  https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz
for c in $(seq 1 22); do
  V=$D/raw/chr$c.vcf.gz
  [[ -e $V.ok ]] || { curl -sS -L -C - --retry 5 -o $V $B/ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz && touch $V.ok; }
  [[ -e $D/panel/chr$c.genotypes.bcf.csi ]] || { echo "$c chr$c" > $D/raw/rename_$c.txt
    bcftools view -m2 -M2 -v snps --threads 4 $V -Ou | bcftools annotate --rename-chrs $D/raw/rename_$c.txt -Ob -o $D/panel/chr$c.genotypes.bcf
    bcftools index $D/panel/chr$c.genotypes.bcf; }
  [[ -e $D/sites/chr$c.vcf.gz.csi ]] || { bcftools view -G -m2 -M2 -v snps -i 'INFO/AF>=0.05' $V -Oz -o $D/sites/chr$c.vcf.gz
    bcftools index $D/sites/chr$c.vcf.gz; }
done
bcftools concat $(for c in $(seq 1 22); do echo $D/sites/chr$c.vcf.gz; done) -Oz -o $D/snps_AF5e2_chr1to22_hg19.vcf.gz
bcftools index $D/snps_AF5e2_chr1to22_hg19.vcf.gz
bcftools --version | head -2 > $D/versions.txt
echo done
