#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Download the REBC-THYR claim 7 files listed in DATA/files.tsv (gdc_manifest.py) from the GDC open-access API into
# DATA/<kind>/<file_id>/<file_name>, md5-check each against the GDC md5sum, and write DATA/manifest.tsv.
# Resumable: a file with a .ok marker is skipped; a mismatch is renamed *.bad_md5 and retried once.
# No file content is read beyond md5sum.
# Usage: sbatch -J rebc_fetch --output=DATA/logs/%x_%j.out fetch.sbatch DATA
set -euo pipefail
DATA=$1
cd "$DATA"
tail -n +2 files.tsv | cut -f1-5 | while IFS=$'\t' read -r kind fid name md5 size; do
  d=$kind/$fid; f=$d/$name
  test -f "$f.ok" && continue
  mkdir -p "$d"
  for try in 1 2; do
    curl -sS --fail --retry 5 --retry-delay 20 -C - -o "$f" "https://api.gdc.cancer.gov/data/$fid" || true
    got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1 || true)
    if [ "$got" = "$md5" ]; then touch "$f.ok"; break; fi
    mv -f "$f" "$f.bad_md5" 2>/dev/null || true
    echo "md5 mismatch try $try: $f" >&2
  done
done
{ printf 'kind\tfile_id\tfile_name\tpath\tbytes\tmd5_gdc\tmd5_ok\n'
  tail -n +2 files.tsv | cut -f1-4 | while IFS=$'\t' read -r kind fid name md5; do
    f=$kind/$fid/$name
    if test -f "$f.ok"; then printf '%s\t%s\t%s\t%s\t%s\t%s\tyes\n' "$kind" "$fid" "$name" "$f" "$(stat -c %s "$f")" "$md5"
    else printf '%s\t%s\t%s\t%s\tNA\t%s\tno\n' "$kind" "$fid" "$name" "$f" "$md5"; fi
  done; } > manifest.tsv
echo "files ok: $(grep -c $'\tyes$' manifest.tsv) / $(($(wc -l < manifest.tsv) - 1))"
