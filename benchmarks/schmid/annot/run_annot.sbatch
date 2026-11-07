#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -J masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 deviation (STATUS claim 5, 2026-11-04): download the two raw BioMart tables once from the Ensembl 113 archive
# (the release current on 2025-03-25, the authors' commit 41dfa38 that set www.ensembl.org "#currently: GRCh38.p14"),
# save query, URL, date and md5, then build the CaSpER and CONICSmat annotation files for each dataset.
# Usage: sbatch --output=OUT/annot_%j.out run_annot.sbatch OUT INPUTS_DIR DS [DS ...]
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 INP=$2; shift 2
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
HOST=https://oct2024.archive.ensembl.org
URL=$HOST/biomart/martservice
R=$OUT/raw; mkdir -p $R
fetch() {   # name, attribute list
  local name=$1; shift
  local attrs=""; for a in "$@"; do attrs="$attrs<Attribute name=\"$a\"/>"; done
  local q='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6" completionStamp="1"><Dataset name="hsapiens_gene_ensembl" interface="default">'"$attrs"'</Dataset></Query>'
  echo "$q" > $R/query_$name.xml
  curl -sS -L --fail -m 3000 --retry 3 --data-urlencode "query=$q" "$URL" -o $R/biomart_$name.tsv
  test "$(tail -n1 $R/biomart_$name.tsv)" = "[success]"
  echo -e "$name\t$URL\t$(date -Is)\t$(md5sum < $R/biomart_$name.tsv | cut -d' ' -f1)\t$(wc -l < $R/biomart_$name.tsv)" >> $R/MANIFEST.tsv
}
echo -e "name\turl\tdate\tmd5\tlines" > $R/MANIFEST.tsv
curl -sS -L -m 60 "$HOST/biomart/martservice?type=version&mart=ENSEMBL_MART_ENSEMBL" > $R/mart_version.txt || true
curl -sS -L -m 60 "$URL?type=registry" > $R/registry.xml || true
fetch casper hgnc_symbol chromosome_name start_position end_position band
fetch conicsmat ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
grep -o 'displayName="Ensembl Genes [0-9]*"' $R/registry.xml | head -1 > $R/release.txt || true
for DS in "$@"; do
  $MM run -n schmid_r41 Rscript $S/annot/build_annotations.R $R $INP/input_$DS/count_matrix.txt $DS $OUT/$DS \
    > $OUT/build_$DS.log 2>&1
done
md5sum $OUT/*/gene_annotation_* > $OUT/md5_annotations.txt
echo ANNOT_DONE
