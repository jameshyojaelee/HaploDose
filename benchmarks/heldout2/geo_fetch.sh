#!/usr/bin/env bash
# Download one GEO supplementary file over HTTPS (resumable) and append a manifest row. GEO publishes no md5, so the
# md5 is computed locally after the size matches the server's Content-Length (md5_source = local_after_size_match).
# Usage: geo_fetch.sh <url> <dest_dir>
set -euo pipefail
u=$1; d=$2; f=$d/$(basename "$u"); mkdir -p "$d"
[[ -e $f.ok ]] && exit 0
n=$(curl -sSIL "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{v=$2} END{print v}')
for i in $(seq 1 200); do
  [[ -e $f && $(stat -c %s "$f") -eq $n ]] && break
  curl -sS -L -C - --retry 10 --retry-delay 30 --speed-limit 1000 --speed-time 300 -o "$f" "$u" || sleep 30
done
[[ $(stat -c %s "$f") -eq $n ]] || { echo "SIZE MISMATCH $f" >&2; exit 1; }
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(basename "$d")/$(basename "$f")" "$u" "$n" "$(md5sum "$f" | cut -d' ' -f1)" \
  local_after_size_match "$(date -Iseconds)" > "$f.ok"
