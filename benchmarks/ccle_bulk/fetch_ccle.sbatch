#!/usr/bin/env bash
#SBATCH --job-name=fetch_ccle
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/%x_%j.out
# DepMap 24Q4 Public (figshare article 27993248) files for the CCLE bulk benchmark (prespecified SOTA claim 3):
# RNA-SeQC gene counts, ABSOLUTE CN segments (TCN, ploidy), relative CN segments, profile/model maps.
# md5 checked against figshare computed_md5; manifest.tsv records url, bytes, md5, time.
set -euo pipefail
D=data/external/depmap_24q4_ccle
mkdir -p $D
printf 'file\turl\tbytes\tmd5\tretrieved\n' > $D/manifest.tsv
while read -r id name md5; do
  u=https://ndownloader.figshare.com/files/$id
  [[ -e $D/$name.ok ]] || { curl -sSL --retry 10 --retry-delay 30 -o $D/$name $u; echo "$md5  $D/$name" | md5sum -c - && touch $D/$name.ok; }
  printf '%s\t%s\t%s\t%s\t%s\n' $name $u $(stat -c %s $D/$name) $md5 $(date -Is) >> $D/manifest.tsv
done <<'L'
51065501 OmicsExpressionRNASeQCGeneCountProfile.csv 9ab6b559a83b0886563cf7f52e5da15e
51065321 OmicsAbsoluteCNSegmentsProfile.csv e5d3a4425a3ce6d77c42ca75cd13b71b
51065333 OmicsCNSegmentsProfile.csv 79fc5cf5c46f215ef085dddaa0e56860
51065339 OmicsDefaultModelProfiles.csv 4d342a03701c0a4c0742fedb56c1de40
51065723 OmicsProfiles.csv 7195088a5fc41eddb3f896c6e984a992
51065297 Model.csv 675210d17675f3517b0ce39a3c274f16
L
