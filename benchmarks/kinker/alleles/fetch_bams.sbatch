#!/usr/bin/env bash
#SBATCH --job-name=kinker_fetch_bams
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020 (SRP279492): submitters' Cell Ranger 1.2 possorted BAMs (hg19, STAR 2.5.1b, CB/UB tags) from ENA,
# 19 pools, 292 GB; md5 checked against ENA submitted_md5 (benchmarks/kinker/alleles/ena_bams.tsv). Resumes.
set -euo pipefail
# 2026-11-02: HTTPS gave 403 pages and SSL EOF truncations, one file at a time (~6 GB/h). Now: ftp://, resume until
# the size matches ENA submitted_bytes, a sub-1 kB partial (an error page) is restarted, 4 files in parallel,
# HNSCC_custom skipped (not one of the Kinker pools used here).
D=data/external/kinker2020_bams
mkdir -p $D
fetch() {  # url md5 bytes dest
  local u=$1 m=$2 n=$3 f=$4 i
  [[ -e $f.ok ]] && return 0
  for i in $(seq 1 2000); do   # FTP stalls often; resume keeps the partial file
    [[ -e $f && $(stat -c %s "$f") -lt 1024 ]] && : > "$f"
    [[ -e $f && $(stat -c %s "$f") -eq $n ]] && break
    curl -sS -C - --retry 10 --retry-delay 30 --speed-limit 1000 --speed-time 300 -o "$f" "ftp://$u" || sleep 30
  done
  echo "$m  $f" | md5sum -c - && touch "$f.ok"
}
export -f fetch
tail -n +2 benchmarks/kinker/alleles/ena_bams.tsv | while IFS=$'\t' read -r run title ftp md5 bytes; do
  [[ $title == HNSCC_custom ]] && continue
  IFS=';' read -r u_bam u_bai <<< "$ftp"; IFS=';' read -r m_bam m_bai <<< "$md5"; IFS=';' read -r n_bam n_bai <<< "$bytes"
  printf '%s\t%s\t%s\t%s\n' "$u_bai" "$m_bai" "$n_bai" "$D/$run.bam.bai" "$u_bam" "$m_bam" "$n_bam" "$D/$run.bam"
done | xargs -P 4 -L 1 bash -c 'fetch "$0" "$1" "$2" "$3"'
ls $D/*.ok | wc -l
