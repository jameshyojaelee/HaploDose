# Score the wrapper output with the Schmid scorer functions (their lib.R, pinned) plus lib_haplodose.R:
# read_wgs_results -> read_haplodose -> combine_range_objects -> Pearson, AUC (ROCR), truncated AUC (pROC).
# This is the per-method path (evaluate_results_permethod.R scores each method on its own bins) without the
# optimal-F1 step, which needs crfsuite (not installed). Also checks read_haplodose against a direct pseudobulk.
# Usage: Rscript score_test.R OUT DATASET
suppressPackageStartupMessages({library(data.table); library(dplyr); library(GenomicRanges); library(ROCR); library(pROC)})
args<-commandArgs(trailingOnly=TRUE)
out<-args[1]; ds<-args[2]
H<-"/gpfs/commons/home/jameslee/Haplodose"
source(file.path(H,"external/benchmark_scrnaseq_cnv_callers/snakemake_pipeline/scripts/lib.R"))
source(file.path(H,"benchmarks/schmid/lib_haplodose.R"))
inp<-file.path(out,paste0("input_",ds))
hd_dir<-file.path(out,"output",ds,"haplodose")
annot<-file.path(inp,"sample_annotation.txt"); refs<-file.path(inp,"ref_groups.txt")
genes_f<-file.path(hd_dir,paste0(ds,"_haplodose_gene_by_clone.tsv"))
cells_f<-file.path(hd_dir,paste0(ds,"_haplodose_cell_clone.tsv"))

wgs<-read_wgs_results(file.path(out,paste0(ds,"_WGS_groundtruth"),"wgs_results_formated.csv"))
hd<-read_haplodose(genes_f,cells_f,annot,refs)
cat("truth bins",length(wgs),"| HaploDose genes",length(hd[[1]]),"| non-reference assigned cells",hd[[2]],"\n")

#Check the pseudobulk against the per-cell reader: the mean over cells must equal the clone-weighted mean
ind<-read_haplodose_individual_cells(genes_f,cells_f,annot,refs)
stopifnot(ncol(ind[[2]])==hd[[2]], max(abs(rowMeans(ind[[2]])-hd[[1]]$haplodose))<1e-9)
#Reference cells must not enter the pseudobulk
a<-fread(annot,header=FALSE); r<-fread(refs)
stopifnot(!any(colnames(ind[[2]]) %in% a$V1[a$V2 %in% r$ref_groups]))

comb<-combine_range_objects(wgs,hd[[1]],method_colname="haplodose")
m<-as.data.frame(elementMetadata(comb))
corr<-evaluate_correlation_colnames(m,"wgs_mean","haplodose",printres=FALSE)
auc<-evaluate_gains_loss_auc(m,"haplodose","wgs_mean")
auct<-evaluate_gains_loss_auc_truncated(m,"haplodose","wgs_mean")
res<-data.frame(method="haplodose",num_bins=length(comb),pearson=corr[1],auc_gains=auc[1],auc_gains_trunc=auct[1],
                aucpr_gains=auc[2],auc_losses=auc[3],auc_losses_trunc=auct[2],aucpr_losses=auc[4],
                n_gain_bins=sum(m$wgs_mean>2.5),n_loss_bins=sum(m$wgs_mean<1.5))
print(res)
fwrite(res,file.path(out,"score_test.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo_score.txt"))
