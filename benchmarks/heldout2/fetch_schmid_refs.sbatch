#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_schmid_refs
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# External reference cells for Schmid 2025 COLO320 and MCF7 (claim 5). Authors' repo commit 319ea34:
#   prepare_COLO320_rep1_withRef.R reads "Adult_healthy_large_intestine.h5ad" (X = raw counts, obs$category ==
#     "Epithelial", 5000 cells sampled with set.seed(1)); the repo has no URL or subsetting code for it. Supplementary
#     Data 1 (MOESM3) names the source: "Healthy human adult from gut cell atlas (Elmentaite et al., 2021)".
#     Fetched: gutcellatlas.org raw-count objects Full_obj_raw_counts_nosoupx_v2.h5ad (1.84 GB, primary candidate)
#     and epi_raw_counts02_v2.h5ad (0.67 GB, epithelial-only alternative). Which one, and the adult/healthy/large-
#     intestine subset, is [Inference] and must be fixed before scoring.
#   prepare_MCF7_dataset_withRef.R reads data/Tabular_sapiens/TS_Mammary.h5ad (layer raw_counts, obs method ==
#     "10X", obs compartment). Fetched: Tabula Sapiens figshare article 14267219 file 34701997 TS_Mammary.h5ad.zip
#     (381 MB; same file in article versions 4 and 5; version 3 had file 28902069). Which version the authors used
#     is unverified.
# Checks: gut atlas against the S3 ETag (multipart md5), figshare against its API computed_md5. Total about 2.9 GB.
set -euo pipefail
G=data/external/gutcellatlas_elmentaite2021; T=data/external/tabula_sapiens_v1_figshare
mkdir -p $G $T
get() {  # url dest bytes
  for i in $(seq 1 100); do
    [[ -e $2 && $(stat -c %s "$2") -eq $3 ]] && return 0
    curl -sS -L -C - --retry 10 --retry-delay 30 --speed-limit 1000 --speed-time 300 -o "$2" "$1" || sleep 30
  done
  return 1
}
hdr='file\turl\tbytes\tmd5\tmd5_source\tretrieved\n'
printf "$hdr" > $G/manifest.tsv
for f in Full_obj_raw_counts_nosoupx_v2.h5ad epi_raw_counts02_v2.h5ad; do
  u=https://cellgeni.cog.sanger.ac.uk/gutcellatlas/$f
  h=$(curl -sSI $u | tr -d '\r'); n=$(awk 'tolower($1)=="content-length:"{print $2}' <<< "$h"); e=$(awk 'tolower($1)=="etag:"{print $2}' <<< "$h")
  get $u $G/$f $n
  read -r m st < <(python3 -I benchmarks/heldout2/s3_etag_check.py $G/$f "$e")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' $f $u $n $m "s3_etag_${e//\"/}:${st#etag:}" "$(date -Iseconds)" >> $G/manifest.tsv
done
u=https://ndownloader.figshare.com/files/34701997; f=TS_Mammary.h5ad.zip; want=d7946751e39eefb4a69d58aade120bf8
get $u $T/$f 380791465
m=$(md5sum $T/$f | cut -d' ' -f1); [[ $m == $want ]] && st=match || st=MISMATCH
printf "$hdr" > $T/manifest.tsv
printf '%s\t%s\t%s\t%s\t%s\t%s\n' $f $u 380791465 $m "figshare_computed_md5:$st" "$(date -Iseconds)" >> $T/manifest.tsv
[[ $st == match ]] && (cd $T && unzip -o -q $f && ls -la)
cat $G/manifest.tsv $T/manifest.tsv
