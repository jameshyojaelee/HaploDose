#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_gse123813
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Yost et al. 2019 BCC/SCC (GEO GSE123813), all 7 supplementary files (~150 MB). Schmid prepare_BCC_Sample06.R uses
# GSE123813_bcc_scRNA_counts.txt.gz and GSE123813_bcc_all_metadata.txt.gz (GEO matrix + published annotation).
# GEO has no md5: size checked against Content-Length, md5 computed locally (benchmarks/heldout2/geo_fetch.sh).
set -euo pipefail
D=data/external/yost2019_gse123813
for f in bcc_all_metadata.txt.gz bcc_scRNA_counts.txt.gz bcc_tcell_metadata.txt.gz bcc_tcr.txt.gz \
         scc_metadata.txt.gz scc_scRNA_counts.txt.gz scc_tcr.txt.gz; do
  bash benchmarks/heldout2/geo_fetch.sh https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_$f $D/suppl
done
{ printf 'file\turl\tbytes\tmd5\tmd5_source\tretrieved\n'; cat $D/suppl/*.ok; } > $D/manifest.tsv
cat $D/manifest.tsv
