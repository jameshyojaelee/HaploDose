#!/usr/bin/env Rscript
# R inferCNV 1.22 (module R/4.4.1 + JAGS/4.3.2, user library 4.4), called by run_infercnv_r.py.
# Settings copied from ~/CRC/src/43_cnv_caller/benchmark_infercnv.R: cutoff 0.1, window 100,
# samples mode, no denoise, no HMM, cells with >= 300 UMI. Intermediate .infercnv_obj files are
# not saved (they reach several GB per section). The final residual matrix (genes x cells,
# fit reference and observation cells) is written as raw float32 for Python.
# args: prep_dir (counts.mtx.gz genes x cells, genes.tsv, cells.tsv, annotations.tsv, gene_order.tsv), out_dir, ncores
set.seed(20261101)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
prep <- args[1]; out <- args[2]; ncores <- as.integer(args[3])
stopifnot(!startsWith(normalizePath(out, mustWork = FALSE), '/gpfs/commons/home/jameslee/CRC'))
.libPaths(c('/gpfs/commons/home/jameslee/R/x86_64-pc-linux-gnu-library/4.4', .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(infercnv)})
dir.create(out, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[['elapsed']]
counts <- as(readMM(gzfile(file.path(prep, 'counts.mtx.gz'))), 'CsparseMatrix')
counts <- as(counts, 'dgCMatrix')
rownames(counts) <- read.delim(file.path(prep, 'genes.tsv'))$gene
colnames(counts) <- read.delim(file.path(prep, 'cells.tsv'))$cell_id
annotation <- read.delim(file.path(prep, 'annotations.tsv'), header = FALSE)
obj <- CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = file.path(prep, 'annotations.tsv'),
                            gene_order_file = file.path(prep, 'gene_order.tsv'),
                            ref_group_names = unique(annotation$V2[grepl('^ref_', annotation$V2)]),
                            min_max_counts_per_cell = c(300, Inf))
rm(counts); invisible(gc())
t1 <- proc.time()[['elapsed']]
obj <- infercnv::run(obj, cutoff = .1, out_dir = out, cluster_by_groups = TRUE,
                     denoise = FALSE, HMM = FALSE, no_plot = TRUE, no_prelim_plot = TRUE,
                     analysis_mode = 'samples', output_format = 'pdf', resume_mode = FALSE,
                     num_threads = ncores, window_length = 100, write_expr_matrix = FALSE,
                     save_rds = FALSE, save_final_rds = FALSE)
t2 <- proc.time()[['elapsed']]
E <- obj@expr.data
writeLines(rownames(E), file.path(out, 'expr_genes.txt'))
writeLines(colnames(E), file.path(out, 'expr_cells.txt'))
con <- file(file.path(out, 'expr_cells_by_genes.f32'), 'wb')
for (lo in seq(1, ncol(E), by = 2000)) {          # cell-major float32, written in blocks
    hi <- min(lo + 1999, ncol(E))
    writeBin(as.vector(E[, lo:hi, drop = FALSE]), con, size = 4)
}
close(con)
write.table(data.frame(method = 'R_infercnv', n_cells = ncol(E), n_genes = nrow(E), load_seconds = t1 - t0,
                       inference_seconds = t2 - t1, HMM = FALSE, analysis_mode = 'samples'),
            file.path(out, 'timing.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
