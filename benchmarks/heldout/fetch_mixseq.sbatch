#!/usr/bin/env bash
#SBATCH --job-name=fetch_mixseq
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout/%x_%j.out
# MIX-seq (McFarland et al. 2020, Figshare 10298696 v3, CC-BY 4.0): control pools only (DMSO, untreated), the
# supplementary tables and README; md5 checked against Figshare's supplied md5 (benchmarks/heldout/mixseq_files.tsv).
# Zips are kept as downloaded; extraction is a separate step into its own directory.
set -euo pipefail
D=data/external/mixseq_mcfarland2020/download
mkdir -p $D
while IFS=$'\t' read -r name url md5 size; do
  f="$D/$name"
  [[ -e $f.ok ]] && continue
  curl -sS -L -C - --retry 10 --retry-delay 30 -o "$f" "$url"
  echo "$md5  $f" | md5sum -c - && touch "$f.ok"
done < benchmarks/heldout/mixseq_files.tsv
ls -la $D
