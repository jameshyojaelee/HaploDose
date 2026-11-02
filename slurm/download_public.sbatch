#!/usr/bin/env bash
#SBATCH --job-name=hd_download
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/download/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/download/%x_%j.err
# Public benchmark data (open GEO supplementary files), one new directory per set, with sha256 and source.
# Usage: download_public.sbatch   (skips sets whose directory already exists)
set -euo pipefail
ROOT=/gpfs/commons/home/jameslee/Haplodose/data/external
mkdir -p "$ROOT"
fetch() {   # name url...
  local name=$1; shift
  local dir=$ROOT/$name
  [[ -e $dir ]] && { echo "exists: $dir"; return; }
  local tmp=$ROOT/.$name.part.$SLURM_JOB_ID
  mkdir -p "$tmp"
  for url in "$@"; do
    wget -q --tries=5 --waitretry=30 -P "$tmp" "$url"
    printf '%s\t%s\t%s\n' "$(basename "$url")" "$url" "$(date -u +%FT%TZ)" >> "$tmp/source.tsv"
  done
  (cd "$tmp" && sha256sum $(cut -f1 source.tsv) > sha256.txt)
  mv -T "$tmp" "$dir"
  echo "done: $dir"; du -sh "$dir"
}
G=https://ftp.ncbi.nlm.nih.gov/geo/series
fetch andor2020_gse142750 $G/GSE142nnn/GSE142750/suppl/GSE142750_RAW.tar $G/GSE142nnn/GSE142750/suppl/filelist.txt
fetch kinker2020_gse157220 $G/GSE157nnn/GSE157220/suppl/GSE157220_CPM_data.txt.gz $G/GSE157nnn/GSE157220/suppl/GSE157220_Pool_composition.xlsx
