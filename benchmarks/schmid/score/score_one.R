# Claim 5 scoring of one dataset in one scoring working dir W (built by run_score_all.sbatch), with the
# HaploDose-patched scorers and the snakemake stand-in. Inputs and params mirror their evaluation rule for the
# dataset (workflow.sm at 319ea34): evaluate_cnv_prediction_wes (BCC06, BCC06post, MM: GATK / WES),
# evaluate_cnv_prediction_cnvkit_COLO320 (COLO320: genomic_format CNVkit, as in their rule) and evaluate_cnv_prediction_wgs
# (MCF7: no params). Gene positions for the evaluation: hg38_gencode_v27_fixtab_20261104.txt (STATUS claim 5).
#   MODE combined  : scripts/evaluate_results_haplodose.R (primary; common bins across all methods incl. HaploDose)
#   MODE permethod : scripts/evaluate_results_permethod_haplodose.R (sensitivity; each method on its own bins;
#                    same genomic_format as the combined rule)
# Usage (from W): Rscript score_one.R SHIM DS MODE
args <- commandArgs(trailingOnly = TRUE)
shim <- args[1]; ds <- args[2]; mode <- args[3]
source(shim)
truth <- switch(ds,
  COLO320 = "data/COLO320_WGS/COLO320_HSR_WGS.aln.sort.rmdup.call.cns",
  MCF7 = "data/MCF7_WGS_groundtruth/wgs_results_formated.csv",
  paste0("data/", ds, "_WES_GATK/", ds, "_tumor_clean.called.seg"))
# COLO320 fix 2026-11-05: their rule evaluate_cnv_prediction_cnvkit_COLO320 sets genomic_format = "CNVkit" (workflow.sm at
# 319ea34); without it the scorer reads the .cns truth as WGS and stops ("object 'loss_wgs' not found").
params <- switch(ds, BCC06 = , BCC06post = , MM = list(genomic_format = "GATK", genomic_name = "WES"),
  COLO320 = list(genomic_format = "CNVkit"), list())
params$haplodose_lib <- "scripts/lib_haplodose.R"
suf <- switch(ds, COLO320 = "_cnvkit", MCF7 = "", "_wes")
o <- paste0("results/output_", ds, "/")
inp <- list(
  annot = paste0("data/input_", ds, "/sample_annotation.txt"),
  ref_groups = paste0("data/input_", ds, "/ref_groups.txt"),
  genomic_ground_truth = truth,
  infercnv_cnv = paste0(o, "infercnv/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"),
  infercnv_expr = paste0(o, "infercnv/infercnv.observations.txt"),
  infercnv_gene_pos = "data/annotations/hg38_gencode_v27_fixtab_20261104.txt",
  casper = paste0(o, "casper/", ds, "_casper_pseudobulk_aggregate.RDS"),
  copykat = paste0(o, "copykat/", ds, "_copykat_CNA_raw_results_gene_by_cell.txt"),
  scevan = paste0(o, "scevan/output/", ds, "_CNAmtxSubclones.RData"),
  scevan_gene_pos = paste0(o, "scevan/output/", ds, "_count_mtx_annot.RData"),
  scevan_clone1 = paste0(o, "scevan/output/", ds, "_subclone1_CN.seg"),
  scevan_annot = paste0(o, "scevan/output/", ds, "_scevan_prediction.txt"),
  numbat_obj = paste0(o, "numbat/gexp_roll_wide.tsv.gz"),
  CONICSmat_filtered_dims = paste0(o, "CONICSmat/", ds, "_filtering_dims.tsv"),
  CONICSmat_chr_pos = "data/annotations/chromosome_arm_positions_grch38.txt",
  CONICSmat_cnv = paste0(o, "CONICSmat/cnv_types.tsv"),
  haplodose_genes = paste0(o, "haplodose/", ds, "_haplodose_gene_by_clone.tsv"),
  haplodose_cells = paste0(o, "haplodose/", ds, "_haplodose_cell_clone.tsv"))
# STATUS 2026-11-04 ~16:45: run_numbat failed on MCF7 (not rerun); score without the Numbat input, as the gate does.
if (ds == "MCF7") inp$numbat_obj <- NULL
miss <- names(inp)[!file.exists(unlist(inp))]
if (length(miss)) stop("missing inputs: ", paste(miss, collapse = ", "))
e <- paste0(o, "evaluation_haplodose/")
dir.create(e, showWarnings = FALSE)
if (mode == "combined") {
  outs <- list(merged_results = "outputs_allmethods_combined", text_filter = "filtering_results",
    plot_filter = "filtering_results", plot_heatmap = "evaluation_cnv_prediction_heatmap",
    text_scaling = "scaling_factors", text_corr = "evaluation_cnv_prediction_corr",
    plot_corr = "evaluation_cnv_prediction_corr", text_auc = "evaluation_cnv_prediction_auc",
    plot_auc = "evaluation_cnv_prediction_auc", text_f1 = "evaluation_cnv_prediction_f1",
    plot_f1 = "evaluation_cnv_prediction_f1", plot_heatmap_binf1 = "evaluation_cnv_prediction_heatmap_binarized_f1")
  ext <- c(merged_results = ".tsv", text_filter = ".tsv", plot_filter = ".pdf", plot_heatmap = ".pdf",
    text_scaling = ".tsv", text_corr = ".tsv", plot_corr = ".pdf", text_auc = ".txt", plot_auc = ".pdf",
    text_f1 = ".txt", plot_f1 = ".pdf", plot_heatmap_binf1 = ".pdf")
  outs <- setNames(lapply(names(outs), function(k) paste0(e, outs[[k]], suf, ext[[k]])), names(outs))
  script <- "scripts/evaluate_results_haplodose.R"
} else if (mode == "permethod") {
  outs <- list(res = paste0(e, "evaluation_cnv_prediction_permethod", suf, ".tsv"))
  script <- "scripts/evaluate_results_permethod_haplodose.R"
} else stop("unknown mode ", mode)
if (any(file.exists(unlist(outs)))) stop("output exists; refusing to overwrite")
snakemake <- make_snakemake(input = inp, output = outs, params = params, wildcards = list(dataset = ds),
                            log = list(paste0("logs/score_", mode, "_", ds, ".log")), rule = paste0("score_", mode))
source(script)
writeLines(capture.output(sessionInfo()), paste0(e, "sessionInfo_", mode, ".txt"))
cat("SCORE_ONE_DONE", ds, mode, "\n")
