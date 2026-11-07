# ------------------------------------------------------------------------------
# HaploDose as one more caller in the Schmid et al. 2025 benchmark pipeline
# (external/benchmark_scrnaseq_cnv_callers, pinned commit in external/README.md).
# Include from their workflow.sm, after copying benchmarks/schmid/run_haplodose.py, lib_haplodose.R and the
# two *_haplodose.R scorers into their snakemake_pipeline/scripts/:
#     include: "haplodose.smk"
# Same inputs as run_copykat. The frozen HaploDose code path is a config value, never the env install:
#     snakemake -s workflow.sm --config haplodose_code=/path/to/frozen/code ...
# CPU by default (haplodose_platform=cpu). Not run under snakemake here (snakemake is not installed); the
# wrapper and readers were tested directly, see benchmarks/schmid/PLAN.md.
# ------------------------------------------------------------------------------

HAPLODOSE_PY = config.get("haplodose_python", "/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python")

rule run_haplodose:
	input:
		matrix = "data/input_{dataset}/count_matrix.txt",
		annot = "data/input_{dataset}/sample_annotation.txt",
		ref_groups = "data/input_{dataset}/ref_groups.txt",
		gene_pos = "data/annotations/hg38_gencode_v27.txt"
	output:
		genes = "results/output_{dataset}/haplodose/{dataset}_haplodose_gene_by_clone.tsv",
		cells = "results/output_{dataset}/haplodose/{dataset}_haplodose_cell_clone.tsv"
	log:
		"logs/haplodose_{dataset}.log"
	benchmark:
		"benchmarks/haplodose_{dataset}_benchmark.txt"
	threads: 8
	resources:
		mem_mb = 32000
	params:
		time = "24:00:00",
		code = config["haplodose_code"],
		platform = config.get("haplodose_platform", "cpu"),
		column = config.get("haplodose_column", "log_ratio_gene"),
		holdout = config.get("haplodose_ref_holdout", 0.3)
	shell:
		"{HAPLODOSE_PY} -I scripts/run_haplodose.py --matrix {input.matrix} --annot {input.annot} "
		"--ref-groups {input.ref_groups} --gene-pos {input.gene_pos} "
		"--out-dir results/output_{wildcards.dataset}/haplodose --dataset {wildcards.dataset} "
		"--code {params.code} --jax-platform {params.platform} --column {params.column} "
		"--ref-holdout {params.holdout} > {log} 2>&1"

# Combined scorer (bins shared by all methods) with HaploDose added: same inputs as their
# evaluate_cnv_prediction_wgs plus the two HaploDose files; script is the patched copy.
# The WES (GATK) and CNVkit (COLO320) variants differ only in genomic_ground_truth and params.genomic_format,
# exactly as their evaluate_cnv_prediction_wes and evaluate_cnv_prediction_cnvkit_COLO320 rules do.

rule evaluate_cnv_prediction_wgs_haplodose:
	input:
		annot = "data/input_{dataset}/sample_annotation.txt",
		ref_groups = "data/input_{dataset}/ref_groups.txt",
		genomic_ground_truth = "data/{dataset}_WGS_groundtruth/wgs_results_formated.csv",
		infercnv_cnv = "results/output_{dataset}/infercnv/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt",
		infercnv_expr = "results/output_{dataset}/infercnv/infercnv.observations.txt",
		infercnv_gene_pos = "data/annotations/hg38_gencode_v27.txt",
		casper = "results/output_{dataset}/casper/{dataset}_casper_pseudobulk_aggregate.RDS",
		copykat = "results/output_{dataset}/copykat/{dataset}_copykat_CNA_raw_results_gene_by_cell.txt",
		scevan = "results/output_{dataset}/scevan/output/{dataset}_CNAmtxSubclones.RData",
		scevan_gene_pos = "results/output_{dataset}/scevan/output/{dataset}_count_mtx_annot.RData",
		scevan_clone1 = "results/output_{dataset}/scevan/output/{dataset}_subclone1_CN.seg",
		scevan_annot = "results/output_{dataset}/scevan/output/{dataset}_scevan_prediction.txt",
		numbat_obj = "results/output_{dataset}/numbat/gexp_roll_wide.tsv.gz",
		CONICSmat_filtered_dims = "results/output_{dataset}/CONICSmat/{dataset}_filtering_dims.tsv",
		CONICSmat_chr_pos = "data/annotations/chromosome_arm_positions_grch38.txt",
		CONICSmat_cnv = "results/output_{dataset}/CONICSmat/cnv_types.tsv",
		haplodose_genes = "results/output_{dataset}/haplodose/{dataset}_haplodose_gene_by_clone.tsv",
		haplodose_cells = "results/output_{dataset}/haplodose/{dataset}_haplodose_cell_clone.tsv"
	output:
		merged_results = "results/output_{dataset}/evaluation_haplodose/outputs_allmethods_combined.tsv",
		text_filter = "results/output_{dataset}/evaluation_haplodose/filtering_results.tsv",
		plot_filter = "results/output_{dataset}/evaluation_haplodose/filtering_results.pdf",
		plot_heatmap = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_heatmap.pdf",
		text_scaling = "results/output_{dataset}/evaluation_haplodose/scaling_factors.tsv",
		text_corr = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_corr.tsv",
		plot_corr = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_corr.pdf",
		text_auc = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_auc.txt",
		plot_auc = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_auc.pdf",
		text_f1 = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_f1.txt",
		plot_f1 = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_f1.pdf",
		plot_heatmap_binf1 = "results/output_{dataset}/evaluation_haplodose/evaluation_cnv_prediction_heatmap_binarized_f1.pdf"
	log:
		"logs/evaluate_prediction_haplodose_{dataset}.log"
	threads: 1
	resources:
		mem_mb = 24000
	params:
		time = "24:00:00",
		haplodose_lib = "scripts/lib_haplodose.R"
	script:
		"scripts/evaluate_results_haplodose.R"
