#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_bian
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Claim 4 secondary, cell level: Bian et al. 2018 CRC scTrio-seq2, GSE97693 / SRP103839 / PRJNA382695 (Smart-seq2-
# like scRNA + single-cell bisulfite (PBAT) DNA from the same cell; DNA CN is to be re-derived from the BS reads).
# Runs (benchmarks/heldout2/bian_runs.tsv): all 1198 RNA-Seq runs (group rna, 614 GB) + the bisulfite runs of the
# 616 cells that also have RNA (group dna_bs, matched on sample title minus the scTrioSeq2Met/Rna prefix, 1.53 TB):
# 1817 runs, 3634 files, 2176 GB. The 675 bisulfite runs without an RNA partner (1.65 TB) are listed in
# bian_runs_dna_unpaired_not_fetched.tsv and are not fetched.
# 2.2 TB may not finish in 72 h; the script resumes, so resubmit it unchanged after a timeout.
# BLIND: download and checksum only.
set -euo pipefail
bash benchmarks/heldout2/ena_fetch.sh benchmarks/heldout2/bian_runs.tsv data/external/bian2018_sctrio 8
