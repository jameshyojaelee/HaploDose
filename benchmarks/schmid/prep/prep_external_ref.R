# Claim 5 input preparation for the external-reference datasets (COLO320, MCF7), in env schmid_r41.
# Mirrors external/benchmark_scrnaseq_cnv_callers/preprocessing_datasets/prepare_COLO320_rep1_withRef.R and
# prepare_MCF7_dataset_withRef.R (commit 319ea34) line for line, with one change: the reference cells are not read
# with zellkonverter::readH5AD (not installed) but from the MatrixMarket subset written by
# extract_reference_counts.py, which holds exactly the cells fixed by STATUS item 5 amendment (f), already filtered
# (COLO320: category Epithelial + 5000 drawn with set.seed(1); MCF7: 10X + epithelial) and in their column order.
# Gene intersection, cbind, annotation and write.table calls are theirs. Also writes barcodes.tsv.gz (the cell-line
# barcodes, i.e. the CellRanger filtered barcodes, for Numbat's pileup).
# Usage: Rscript prep_external_ref.R DATASET CELLRANGER_FILTERED_DIR REF_SUBSET_DIR OUT_INPUT_DIR
library(Matrix)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]; cr_dir <- args[2]; ref_dir <- args[3]; out_dir <- args[4]
stopifnot(dataset %in% c("COLO320", "MCF7"), !dir.exists(out_dir))

# Their code: CellRanger filtered matrix with gene symbols as rownames
counts <- readMM(file.path(cr_dir, "matrix.mtx.gz"))
genenames <- fread(file.path(cr_dir, "features.tsv.gz"), header = FALSE)
cells <- fread(file.path(cr_dir, "barcodes.tsv.gz"), header = FALSE)
colnames(counts) <- cells$V1
rownames(counts) <- genenames$V2
print(dim(counts))

# Reference subset (replaces readH5AD + the category/method/compartment filters + sample())
raw_counts_ref <- readMM(file.path(ref_dir, "matrix.mtx"))
ref_genes <- fread(file.path(ref_dir, "genes.tsv"), header = FALSE)$V1
ref_cells <- fread(file.path(ref_dir, "cells.tsv"), header = FALSE)
rownames(raw_counts_ref) <- ref_genes
colnames(raw_counts_ref) <- ref_cells$V1

intersect_genes <- intersect(rownames(raw_counts_ref), rownames(counts))
print(paste("Intersect genes:", length(intersect_genes)))
raw_counts_ref <- raw_counts_ref[intersect_genes, ]
counts <- counts[intersect_genes, ]

counts_combined <- cbind(raw_counts_ref, counts)
if (dataset == "COLO320") {
  stopifnot(ncol(raw_counts_ref) == 5000, all(ref_cells$V2 == "Epithelial"))
  meta_data_combined <- data.frame(barcode = colnames(counts_combined),
                                   sample = c(rep("Epithelial", 5000), rep("COLO320", ncol(counts))))
  df_refs <- data.frame(ref_groups = setdiff(unique(meta_data_combined$sample), "COLO320"))
} else {
  meta_data_combined <- data.frame(barcode = colnames(counts_combined),
                                   sample = c(as.character(ref_cells$V2), rep("MCF7", ncol(counts))))
  df_refs <- data.frame(ref_groups = as.character(unique(ref_cells$V2)))
}

dir.create(out_dir, recursive = TRUE)
write.table(as.matrix(counts_combined), file = file.path(out_dir, "count_matrix.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(meta_data_combined, file = file.path(out_dir, "sample_annotation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(df_refs, file = file.path(out_dir, "ref_groups.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
file.copy(file.path(cr_dir, "barcodes.tsv.gz"), file.path(out_dir, "barcodes.tsv.gz"))
print(table(meta_data_combined$sample))
print(sum(duplicated(colnames(counts_combined))))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "..", paste0("sessionInfo_", dataset, ".txt")))
