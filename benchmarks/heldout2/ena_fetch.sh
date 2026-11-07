#!/usr/bin/env bash
# Resumable ENA fastq download with md5 check, shared by the heldout2 fetch_*.sbatch jobs.
# Usage: ena_fetch.sh <run_table.tsv> <dest_dir> <parallel>
#   run_table.tsv: header row, then columns run, group, sample_title, read_count, base_count,
#                  fastq_ftp (;-separated), fastq_md5 (;-separated), fastq_bytes (;-separated), as written from the
#                  ENA filereport API (result=read_run).
# Layout: <dest_dir>/fastq/<group>/<run>/<file>.fastq.gz; a <file>.ok marker holds that file's manifest row.
# At the end <dest_dir>/manifest.tsv is rebuilt from all .ok markers (file, url, bytes, md5, md5_source, retrieved).
# Blinding: this script only downloads and checksums; it never reads the reads' content.
# Template: benchmarks/kinker/alleles/fetch_bams.sbatch (ftp://, curl -C - resume until the size matches ENA bytes,
# sub-1 kB partials restarted as error pages).
set -euo pipefail
table=$1; dest=$2; par=${3:-4}
mkdir -p "$dest/fastq"

fetch() {  # url md5 bytes dest
  local u=$1 m=$2 n=$3 f=$4 i got
  [[ -e $f.ok ]] && return 0
  mkdir -p "$(dirname "$f")"
  for i in $(seq 1 500); do
    [[ -e $f && $(stat -c %s "$f") -lt 1024 && $n -ge 1024 ]] && : > "$f"
    [[ -e $f && $(stat -c %s "$f") -eq $n ]] && break
    [[ -e $f && $(stat -c %s "$f") -gt $n ]] && : > "$f"
    curl -sS -C - --retry 10 --retry-delay 30 --speed-limit 1000 --speed-time 300 -o "$f" "ftp://$u" || sleep 30
  done
  got=$(md5sum "$f" | cut -d' ' -f1)
  if [[ $got == "$m" ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${f#$DEST/}" "ftp://$u" "$(stat -c %s "$f")" "$got" "ena_fastq_md5:match" \
      "$(date -Iseconds)" > "$f.ok"
  else
    echo "MD5 MISMATCH $f got $got want $m" >&2
    mv "$f" "$f.bad_md5"   # keep for inspection; a rerun downloads afresh
    return 0
  fi
}
export -f fetch
export DEST=$dest

tail -n +2 "$table" | while IFS=$'\t' read -r run group title rc bc ftp md5 bytes; do
  IFS=';' read -ra U <<< "$ftp"; IFS=';' read -ra M <<< "$md5"; IFS=';' read -ra N <<< "$bytes"
  for k in "${!U[@]}"; do
    printf '%s\t%s\t%s\t%s\n' "${U[$k]}" "${M[$k]}" "${N[$k]}" "$dest/fastq/$group/$run/$(basename "${U[$k]}")"
  done
done | xargs -P "$par" -L 1 bash -c 'fetch "$0" "$1" "$2" "$3"'

want=$(tail -n +2 "$table" | awk -F'\t' '{n+=split($6,a,";")} END{print n}')
have=$(find "$dest/fastq" -name '*.ok' | wc -l)
{ printf 'file\turl\tbytes\tmd5\tmd5_source\tretrieved\n'; find "$dest/fastq" -name '*.ok' | sort | xargs cat; } \
  > "$dest/manifest.tsv"
bad=$(find "$dest/fastq" -name '*.bad_md5' | wc -l)
echo "files verified: $have / $want; md5 mismatches this run: $bad"
[[ $have -eq $want ]]
