#!/usr/bin/env bash
#SBATCH --job-name=h2_dntr_manifest
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Runs after the Schmid fetch (afterany): writes the GSE144296 (= SRP245376) manifest as the DNTRseq_ subset of the
# Schmid manifest, and reports verified vs expected file counts (expected: 2 per run in dntrseq_gse144296_runs.tsv).
set -euo pipefail
D=data/external/schmid2025_raw; O=data/external/dntrseq_gse144296/manifest.tsv
{ head -1 $D/manifest.tsv; grep -P "^fastq/DNTRseq_" $D/manifest.tsv; } > $O
want=$(tail -n +2 benchmarks/heldout2/dntrseq_gse144296_runs.tsv | awk -F'\t' '{n+=split($6,a,";")} END{print n}')
echo "GSE144296 verified files: $(($(wc -l < $O)-1)) / $want"
