#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5: build CaSpER/CONICSmat annotations for one dataset from raw Ensembl 113 BioMart tables already downloaded by
# run_annot.sbatch (no new download). Usage: sbatch --output=OUT/annot_%j.out run_annot_reuse.sbatch OUT RAW_DIR INPUT_DIR DS
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 RAW=$2 INP=$3 DS=$4
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
test ! -e $OUT/$DS
cp $RAW/MANIFEST.tsv $OUT/raw_MANIFEST_reused.tsv; echo "$RAW" > $OUT/raw_dir_reused.txt
(cd $RAW && md5sum biomart_casper.tsv biomart_conicsmat.tsv) > $OUT/raw_md5_check.txt
/gpfs/commons/home/jameslee/.local/bin/micromamba run -n schmid_r41 Rscript $S/annot/build_annotations.R \
  $RAW $INP/count_matrix.txt $DS $OUT/$DS > $OUT/build_$DS.log 2>&1
md5sum $OUT/$DS/gene_annotation_* > $OUT/md5_annotations.txt
echo ANNOT_DONE
