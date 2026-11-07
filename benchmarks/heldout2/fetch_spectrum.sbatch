#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_spectrum
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Claim 4 primary (docs/STATUS.md, prespecified SOTA claims): MSK SPECTRUM ovarian 10x scRNA reads, GSE180661 /
# SRP329496 / PRJNA749041. Runs = CD45-negative (tumour) fraction for the 32 patients that also have DLP+ scDNA on
# Synapse syn66366960 (31 patients CD45N + SPECTRUM-OV-003, which has only UNSORTED runs): 224 runs, 448 fastq files,
# 839.4 GB (benchmarks/heldout2/spectrum_runs.tsv; patient map spectrum_patients.tsv). CD45P runs are not fetched.
# Also the GEO processed matrix (GSE180661_matrix.h5, 30 GB) and cell table (GSE180661_GEO_cells.tsv.gz).
# BLIND: download and checksum only. Do not open for analysis until v-next is frozen (STATUS claim 4).
# Resumable: resubmit the same script after a timeout; verified files are skipped.
set -euo pipefail
D=data/external/spectrum_ov_scrna
bash benchmarks/heldout2/ena_fetch.sh benchmarks/heldout2/spectrum_runs.tsv $D 6 || echo "ENA part incomplete"
for f in GSE180661_GEO_cells.tsv.gz GSE180661_matrix.h5; do
  bash benchmarks/heldout2/geo_fetch.sh https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180661/suppl/$f $D/geo
done
cat $D/geo/*.ok >> $D/manifest.tsv
echo done
