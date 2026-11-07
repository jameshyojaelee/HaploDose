#!/usr/bin/env bash
#SBATCH --job-name=fetch_gtex
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/%x_%j.out
# GTEx v10 RNA-SeQC 2.4.2 gene read counts and sample attributes: tissue-matched normal panel for the CCLE bulk
# benchmark (SOTA claim 3; panel chosen by the user 2026-11-02). md5 checked against the GCS x-goog-hash.
set -euo pipefail
D=data/external/gtex_v10
mkdir -p $D
printf 'file\turl\tbytes\tmd5_base64\tretrieved\n' > $D/manifest.tsv
B=https://storage.googleapis.com/adult-gtex
while read -r path md5; do
  f=$D/$(basename $path)
  [[ -e $f.ok ]] || { curl -sSL --retry 10 --retry-delay 30 -o $f $B/$path
    [[ $(openssl md5 -binary $f | base64) == "$md5" ]] && touch $f.ok || { echo "md5 mismatch $f"; exit 1; }; }
  printf '%s\t%s\t%s\t%s\t%s\n' $(basename $path) $B/$path $(stat -c %s $f) $md5 $(date -Is) >> $D/manifest.tsv
done <<'L'
bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz +EWiWJNdzeBD9VR/GWlBtw==
annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt TEXiXD/F7JKpfaDUvLoSBw==
L
