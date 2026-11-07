# Claim 5: MM199 inputs for the authors' pipeline. The evaluated code of their
# preprocessing_datasets/prepare_MM_sciCNV_scran.Rmd (commit 319ea34), copied in order, with these recorded deviations:
#  - paths: the CellRanger directory and the output directory come from the command line;
#  - skipped: every plot and every use of input_MM_old (old_annot), which is not public and only feeds plots/tables;
#    the eval=FALSE GEO chunk; the last section (compares their own copyKat/SCEVAN outputs, needs pipeline results);
#  - added: set.seed(20261103) before any code (their Rmd sets none; Seurat/scran calls below use their own defaults);
#  - added: a per-cluster marker table (MZB1, CD38, IGKC, B2M) written next to the inputs, so the hard-coded
#    control clusters c(4, 6, 8) can be checked against our CellRanger 8.0.1 / Seurat 5.0.1 clustering.
#  - 2026-11-04: optional CONTROL_CLUSTERS (default "4,6,8", theirs) for the STATUS claim 5 MM199 control-cluster rule;
#    per-cell clusters are written to OUT/../MM_cell_clusters.tsv.
# Usage: Rscript prepare_MM.R CR_FILTERED_DIR OUT_INPUT_DIR [CONTROL_CLUSTERS]
suppressMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(scran)
})
args <- commandArgs(trailingOnly = TRUE)
cr <- args[1]; out <- args[2]
ctrl <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(4,6,8)
cat("control clusters:", ctrl, "\n")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261103)

# Rmd lines 89-103
counts<-readMM(file.path(cr, "matrix.mtx.gz"))
genenames<-read.table(file.path(cr, "features.tsv.gz"))
cellnames<-read.table(file.path(cr, "barcodes.tsv.gz"))

colnames(counts)<-cellnames$V1
rownames(counts)<-genenames$V2

#Combine duplicate gene names to one row
counts<-apply(counts, 2, tapply, as.factor(genenames$V2),sum, na.rm=T)

#Create a Seurat object
seurat_obj <- CreateSeuratObject(counts=counts,min.cells=3,min.features=200)
print(seurat_obj)
rm(counts); gc()

# Rmd lines 126, 140
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 3000 & percent.mt < 5)
print(seurat_obj)

# Rmd lines 151-156
sce <- SingleCellExperiment(list(counts=as.matrix(seurat_obj@assays$RNA@layers$counts)))
clusters <- quickCluster(sce, min.size=100)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce,log=FALSE)
log_counts<-log(normcounts(sce)+1)
seurat_obj[["RNA"]]$data<-as(log_counts,"dgCMatrix")
rm(sce, log_counts); gc()

# Rmd lines 178, 194-195, 202-203, 212-213, 218
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj),
                     verbose=FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Rmd lines 249-250
seurat_obj$large_clusters<-ifelse(seurat_obj$RNA_snn_res.0.5 %in% ctrl,
                                  "control","MM")
print(table(seurat_obj$RNA_snn_res.0.5, seurat_obj$large_clusters))

# Added check (not in their Rmd): marker expression per cluster, log scran-normalised data
markers <- intersect(c("MZB1", "CD38", "IGKC", "B2M", "SDC1", "PTPRC"), rownames(seurat_obj))
dat <- GetAssayData(seurat_obj, layer = "data")[markers, , drop = FALSE]
cl <- seurat_obj$RNA_snn_res.0.5
chk <- do.call(rbind, lapply(levels(cl), function(k) {
  idx <- which(cl == k)
  data.frame(cluster = k, n = length(idx), label = ifelse(k %in% ctrl, "control", "MM"),
             t(sapply(markers, function(g) mean(dat[g, idx] > 0))), check.names = FALSE)
}))
colnames(chk)[-(1:3)] <- paste0("frac_pos_", markers)
write.table(chk, file.path(out, "..", "MM_cluster_marker_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(chk)

# Rmd lines 259-268
cluster_markers <- FindMarkers(seurat_obj, ident.1 = ctrl,
                                 min.pct = 0.25,
                               test.use="wilcox")
print(head(cluster_markers))
de_genes<-cluster_markers[cluster_markers$p_val_adj<0.05,]
print(nrow(de_genes))

# Rmd lines 298-313
meta_data<-seurat_obj@meta.data
write.table(data.frame(barcode=rownames(meta_data), cluster=meta_data$RNA_snn_res.0.5),
            file.path(out, "..", "MM_cell_clusters.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
sample_annot<-data.frame(barcode=rownames(meta_data),
                               sample=meta_data$large_clusters)

write.table(sample_annot, file=file.path(out, "sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(sample_annot$sample),"MM"))
write.table(df_refs, file=file.path(out, "ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

write.table(as.matrix(seurat_obj[["RNA"]]$counts),
            file=file.path(out, "count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "..", "sessionInfo_MM.txt"))
cat("PREPARE_MM_DONE\n")
