#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -J masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5: make the authors' pipeline inputs (count_matrix.txt, sample_annotation.txt, ref_groups.txt,
# barcodes.tsv.gz) for one dataset in OUT/input_<DS>/. Their preprocessing scripts run from a working copy; the pinned
# clone is not modified. Deviations (recorded in OUT/deviations_<DS>.txt):
#  BCC: their prepare_BCC_Sample06.R with only `library(zellkonverter)` removed (not installed; unused in the active
#       code). barcodes.tsv.gz = the dataset's cell IDs from sample_annotation.txt (the authors do not say how their
#       file was made); overlap with our CellRanger filtered barcodes is reported.
#  COLO320/MCF7: prep_external_ref.R (their code with the reference read from the amendment-(f) cell subset).
#  MM: prepare_MM.R (2026-11-04), see that file's header.
# Usage: sbatch --output=OUT/logs/prep_<DS>_%j.out prep_inputs.sbatch OUT DS CR_OUTS_DIR
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2 CR=$3
H=/gpfs/commons/home/jameslee/Haplodose
SRC=$H/external/benchmark_scrnaseq_cnv_callers
P=$S/prep
SEL=$H/results/schmid_reference_cells_20261102_2119
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
test "$(git -C $SRC rev-parse HEAD)" = 319ea34ddfc46bdd4c0cb1ed534f4ece0f5c44ca
mkdir -p $OUT/logs
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ndataset\t$DS\nstart\t$(date -Is)" > $OUT/job_prep_$DS.tsv
case $DS in
  BCC)
    W=$OUT/work_BCC; mkdir -p $W/data/BCC_processed_matrices/input_BCC06 $W/data/BCC_processed_matrices/input_BCC06post
    for f in GSE123813_bcc_all_metadata.txt.gz GSE123813_bcc_scRNA_counts.txt.gz; do
      test -e $H/data/external/yost2019_gse123813/suppl/$f.ok
      ln -s $H/data/external/yost2019_gse123813/suppl/$f $W/data/BCC_processed_matrices/$f
    done
    grep -v '^library(zellkonverter)$' $SRC/preprocessing_datasets/prepare_BCC_Sample06.R > $W/prepare_BCC_Sample06.R
    diff $SRC/preprocessing_datasets/prepare_BCC_Sample06.R $W/prepare_BCC_Sample06.R > $OUT/deviations_BCC.diff || true
    (cd $W && /usr/bin/time -v $MM run -n schmid_r41 Rscript prepare_BCC_Sample06.R) > $OUT/logs/prepare_BCC.log 2>&1
    $MM run -n schmid_r41 Rscript -e 'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[1])' $OUT/sessionInfo_BCC.txt
    for d in BCC06 BCC06post; do
      mv $W/data/BCC_processed_matrices/input_$d $OUT/input_$d
      cut -f1 $OUT/input_$d/sample_annotation.txt | gzip -n > $OUT/input_$d/barcodes.tsv.gz
      n_ann=$(zcat $OUT/input_$d/barcodes.tsv.gz | wc -l)
      n_cr=$(zcat $CR/$d/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l)
      n_both=$(comm -12 <(zcat $OUT/input_$d/barcodes.tsv.gz | sort) <(zcat $CR/$d/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sort) | wc -l)
      echo -e "$d\tannotation_cells\t$n_ann\tcellranger_filtered\t$n_cr\toverlap\t$n_both" >> $OUT/barcode_overlap.tsv
      cut -f2 $OUT/input_$d/sample_annotation.txt | sort | uniq -c > $OUT/input_$d.groups.txt
    done ;;
  COLO320|MCF7)
    if [ $DS = COLO320 ]; then H5=$H/data/external/gutcellatlas_elmentaite2021/Full_obj_raw_counts_nosoupx_v2.h5ad; LAYER=X
    else H5=$H/data/external/tabula_sapiens_v1_figshare/TS_Mammary.h5ad; LAYER=raw_counts; fi
    /usr/bin/time -v $MM run -n haplodose python -I $P/extract_reference_counts.py $H5 $LAYER \
      $SEL/${DS}_reference_cells.tsv $OUT/ref_subset_$DS > $OUT/logs/extract_$DS.log 2>&1
    /usr/bin/time -v $MM run -n schmid_r41 Rscript $P/prep_external_ref.R $DS $CR/$DS/outs/filtered_feature_bc_matrix \
      $OUT/ref_subset_$DS $OUT/input_$DS > $OUT/logs/prepare_$DS.log 2>&1
    $MM run -n haplodose pip freeze > $OUT/pip_freeze_haplodose.txt
    cut -f2 $OUT/input_$DS/sample_annotation.txt | sort | uniq -c > $OUT/input_$DS.groups.txt ;;
  MM)
    # prepare_MM.R: the evaluated code of their prepare_MM_sciCNV_scran.Rmd minus plots and input_MM_old (see its header).
    # barcodes.tsv.gz = CellRanger filtered barcodes, as for COLO320/MCF7 (the other CellRanger-based inputs).
    /usr/bin/time -v $MM run -n schmid_r41 Rscript $P/prepare_MM.R $CR/MM/outs/filtered_feature_bc_matrix $OUT/input_MM ${CONTROL_CLUSTERS:-4,6,8} \
      > $OUT/logs/prepare_MM.log 2>&1
    cp $CR/MM/outs/filtered_feature_bc_matrix/barcodes.tsv.gz $OUT/input_MM/barcodes.tsv.gz
    diff <(sed -n 87,316p $SRC/preprocessing_datasets/prepare_MM_sciCNV_scran.Rmd) $P/prepare_MM.R > $OUT/deviations_MM.diff || true
    cut -f2 $OUT/input_MM/sample_annotation.txt | sort | uniq -c > $OUT/input_MM.groups.txt ;;
  *) echo "unknown dataset $DS"; exit 2 ;;
esac
md5sum $OUT/input_*/* > $OUT/md5_inputs_$DS.txt 2>/dev/null || true
echo -e "end\t$(date -Is)" >> $OUT/job_prep_$DS.tsv
echo PREP_DONE
