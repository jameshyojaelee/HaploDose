#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -J masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 5 gate (STATUS item 5, amendment d): rerun the six comparators in the authors' pipeline (pinned commit,
# env schmid_r41) on ONE primary dataset and score them with their unpatched scorer. HaploDose is not run here.
# Submitted 2026-11-03 on team-lead's go (BCC06, BCC06post, COLO320; MCF7 as a line-overlap comparator run).
# Their rules run locally inside this one job (no cluster submission, no arrays): snakemake --cores 8 with a
# 150 GB memory budget, so at most two 64 GB CaSpER-size rules or four 32 GB rules run at once.
# Usage: sbatch --output=OUT/gate_%j.out benchmarks/schmid/gate/run_gate.sbatch OUT DATASET INPUT_DIR TRUTH BAM_DIR
# Optional env (2026-11-04 fixes, STATUS claim 5): ANNOT_DIR = dir with <DS>/gene_annotation_{casper,conicsmat}_<DS>.tsv
#   built from Ensembl 113 archive tables (placed as their rule outputs); NO_NUMBAT=1 = skip Numbat and score with a
#   no-Numbat copy of their evaluation rule (provisional rows). RETICULATE_PYTHON points inferCNV's leiden step at the
#   env python. The job exits non-zero if snakemake fails; it writes OUT/gate_comparison.tsv when scoring ran.
#   DATASET   COLO320 | BCC06 | BCC06post | MM (the authors' dataset names); MCF7 (line overlap, WGS rule)
#   INPUT_DIR their data/input_<DATASET>/ made by their preprocessing script (count_matrix.txt, sample_annotation.txt,
#             ref_groups.txt, barcodes.tsv.gz)
#   TRUTH     the truth file for DATASET in data/external/schmid2025_truth (linked, never read by this script)
#   BAM_DIR   CellRanger outs/ with possorted_genome_bam.bam(.bai) (CaSpER, Numbat)
set -euo pipefail
# Scripts this job calls come from CODE_DIR (a copy made at submission by benchmarks/schmid/submit_frozen.sh), else the live tree.
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2 INP=$3 TRUTH=$4 BAMD=$5
H=/gpfs/commons/home/jameslee/Haplodose
SRC=$H/external/benchmark_scrnaseq_cnv_callers
E=$H/external/schmid_env
A=$OUT/pipeline/data/annotations
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
# inferCNV leiden via reticulate: a job-local symlink to the env python, so reticulate does not treat it as a conda
# python and look for ~/.local/bin/conda (taken from the env's conda-meta/history; STATUS claim 5, 2026-11-04)
mkdir -p $OUT/pyshim/bin; ln -sfn $MAMBA_ROOT_PREFIX/envs/schmid_r41/bin/python3.10 $OUT/pyshim/bin/python
export RETICULATE_PYTHON=$OUT/pyshim/bin/python
ANNOT_DIR=${ANNOT_DIR:-} NO_NUMBAT=${NO_NUMBAT:-0}
test "$(git -C $SRC rev-parse HEAD)" = 319ea34ddfc46bdd4c0cb1ed534f4ece0f5c44ca
test ! -e $OUT/pipeline
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ndataset\t$DS\nstart\t$(date -Is)" > $OUT/job.tsv

# Working copy of their pipeline (the pinned clone stays untouched)
mkdir -p $OUT/pipeline/data $OUT/pipeline/tools $A
cp -r $SRC/snakemake_pipeline/workflow.sm $SRC/snakemake_pipeline/scripts $SRC/snakemake_pipeline/config $OUT/pipeline/
# CaSpER under Seurat 5.0.1: their casper.R reads @assays$RNA@data (Seurat v3/v4 Assay). One line after library(Seurat)
# makes CreateSeuratObject build a v3 Assay (approved 2026-11-04). Original kept as casper.R.orig; diff in OUT.
C=$OUT/pipeline/scripts/casper.R; cp $C $C.orig
test "$(grep -c '^library(Seurat)$' $C)" = 1
sed -i 's/^library(Seurat)$/library(Seurat)\noptions(Seurat.object.assay.version = "v3")/' $C
diff $C.orig $C > $OUT/casper_patch.diff || true
test "$(grep -c '^options(Seurat.object.assay.version = "v3")$' $C)" = 1
# Evaluation rules read gene positions with data.table::fread, which stops at line 82 of hg38_gencode_v27.txt (one trailing
# tab, TNFRSF18). Their evaluation rules (infercnv_gene_pos) point at a copy with that tab stripped; the method rules
# keep the original (read.table reads both identically). STATUS claim 5, 2026-11-04. Diff in OUT/workflow_fixtab.diff.
W=$OUT/pipeline/workflow.sm; cp $W $W.orig
sed -i 's#infercnv_gene_pos = "data/annotations/hg38_gencode_v27.txt"#infercnv_gene_pos = "data/annotations/hg38_gencode_v27_fixtab_20261104.txt"#' $W
diff $W.orig $W > $OUT/workflow_fixtab.diff || true
test "$(grep -c 'hg38_gencode_v27_fixtab_20261104.txt' $W)" = 14
D=$OUT/pipeline/data/input_$DS; mkdir -p $D            # links only; INPUT_DIR and BAM_DIR are not modified
for f in count_matrix.txt sample_annotation.txt ref_groups.txt barcodes.tsv.gz; do ln -s $INP/$f $D/$f; done
ln -s $BAMD/possorted_genome_bam.bam $D/possorted_genome_bam.bam
ln -s $BAMD/possorted_genome_bam.bam.bai $D/possorted_genome_bam.bam.bai
ln -s $E/tools/BAFExtract $OUT/pipeline/tools/BAFExtract
ln -s $E/tools/Eagle_v2.4.1 $OUT/pipeline/tools/Eagle_v2.4.1
ln -s $H/external/annotations/hg38_gencode_v27.txt $A/hg38_gencode_v27.txt
ln -s $H/external/annotations/hg38_gencode_v27_fixtab_20261104.txt $A/hg38_gencode_v27_fixtab_20261104.txt
# Gate reference inputs: pinned dir from fetch_refs.sbatch + fetch_refs2.sbatch (md5 in its MANIFEST.tsv);
# a missing one fails early
RF=$E/ref_20261102_2215
for f in chromosome_arm_positions_grch38.txt genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz genetic_map_hg38_withX.txt.gz \
         1000G_hg38 genome_input_BAFExtract_hg38/hg38_genome_list genome_input_BAFExtract_hg38/hg38/1.bin; do
  test -e $RF/$f || { echo "MISSING reference input $RF/$f"; exit 2; }
done
for f in chromosome_arm_positions_grch38.txt genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz genetic_map_hg38_withX.txt.gz \
         1000G_hg38 genome_input_BAFExtract_hg38; do
  ln -s $RF/$f $A/$f
done
case $DS in
  COLO320) mkdir -p $OUT/pipeline/data/COLO320_WGS; ln -s $TRUTH $OUT/pipeline/data/COLO320_WGS/COLO320_HSR_WGS.aln.sort.rmdup.call.cns
           TARGET=results/output_COLO320/evaluation/evaluation_cnv_prediction_auc_cnvkit.txt; RULE=evaluate_cnv_prediction_cnvkit_COLO320; SUF=_cnvkit ;;
  BCC06|BCC06post|MM) mkdir -p $OUT/pipeline/data/${DS}_WES_GATK; ln -s $TRUTH $OUT/pipeline/data/${DS}_WES_GATK/${DS}_tumor_clean.called.seg
           TARGET=results/output_$DS/evaluation/evaluation_cnv_prediction_auc_wes.txt; RULE=evaluate_cnv_prediction_wes; SUF=_wes ;;
  MCF7) mkdir -p $OUT/pipeline/data/MCF7_WGS_groundtruth; ln -s $TRUTH $OUT/pipeline/data/MCF7_WGS_groundtruth/wgs_results_formated.csv
           TARGET=results/output_MCF7/evaluation/evaluation_cnv_prediction_auc.txt; RULE=evaluate_cnv_prediction_wgs; SUF= ;;
  *) echo "unknown dataset $DS"; exit 2 ;;
esac
for f in count_matrix.txt sample_annotation.txt ref_groups.txt barcodes.tsv.gz possorted_genome_bam.bam possorted_genome_bam.bam.bai; do
  test -e $D/$f || { echo "MISSING input $f"; exit 2; }
done

if [ -n "$ANNOT_DIR" ]; then
  mkdir -p $OUT/pipeline/results/output_$DS/casper $OUT/pipeline/results/output_$DS/CONICSmat
  cp $ANNOT_DIR/$DS/gene_annotation_casper_$DS.tsv $OUT/pipeline/results/output_$DS/casper/
  cp $ANNOT_DIR/$DS/gene_annotation_conicsmat_$DS.tsv $OUT/pipeline/results/output_$DS/CONICSmat/
  echo -e "annot_dir\t$ANNOT_DIR" >> $OUT/job.tsv
fi
EVD=evaluation
if [ "$NO_NUMBAT" = 1 ]; then
  $MM run -n schmid_r41 python -I $S/gate/make_nonumbat_rule.py $OUT/pipeline/workflow.sm $RULE > $OUT/nonumbat_rule.txt
  EVD=evaluation_nonumbat; TARGET=${TARGET/\/evaluation\//\/evaluation_nonumbat\/}
  echo -e "no_numbat\t1" >> $OUT/job.tsv
fi
echo -e "target\t$TARGET" >> $OUT/job.tsv
cd $OUT/pipeline
rc=0
$MM run -n schmid_r41 snakemake -s workflow.sm --cores 8 --resources mem_mb=150000 -k --rerun-incomplete \
  --rerun-triggers mtime --printshellcmds $TARGET > $OUT/snakemake.log 2>&1 || rc=$?
$MM env export -n schmid_r41 > $OUT/env_export.yml
if [ $rc -ne 0 ]; then echo -e "snakemake\tnonzero_exit_$rc" >> $OUT/job.tsv; echo GATE_SNAKEMAKE_FAILED; exit 1; fi
$MM run -n schmid_r41 python -I $S/gate/compare_published.py $OUT/pipeline/results/output_$DS/$EVD $DS "$SUF" \
  $S/gate/published_values.tsv $OUT > $OUT/compare_published.log 2>&1
echo -e "end\t$(date -Is)" >> $OUT/job.tsv
echo GATE_RUN_DONE
