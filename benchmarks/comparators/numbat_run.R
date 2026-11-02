#!/usr/bin/env Rscript
# Numbat on one specimen, called by run_numbat.py. Pattern from ~/CRC/src/42_meeting_followups/T2_run.R.
# args: prep_dir (query.mtx.gz, ref.mtx.gz, genes.tsv, query_cells.tsv, ref_cells.tsv, alleles.tsv.gz), out_dir, ncores,
#       [genome: hg38 (default) or hg19]
set.seed(20261101)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(3, 4))
prep <- args[1]; out <- args[2]; ncores <- as.integer(args[3])
genome <- if (length(args) == 4) args[4] else 'hg38'
stopifnot(genome %in% c('hg38', 'hg19'))
stopifnot(!startsWith(normalizePath(out, mustWork = FALSE), '/gpfs/commons/home/jameslee/CRC'))
suppressPackageStartupMessages({library(numbat); library(Matrix); library(data.table)})
genes <- fread(file.path(prep, 'genes.tsv'))$gene
read_mat <- function(name, cells) {
    m <- as(readMM(gzfile(file.path(prep, name))), 'CsparseMatrix')
    m <- as(m, 'dgCMatrix')
    dimnames(m) <- list(genes, cells)
    m
}
qcells <- fread(file.path(prep, 'query_cells.tsv'))
rcells <- fread(file.path(prep, 'ref_cells.tsv'))
count_mat <- read_mat('query.mtx.gz', qcells$cell)
ref_mat <- read_mat('ref.mtx.gz', rcells$cell)
reference <- aggregate_counts(ref_mat, data.frame(cell = rcells$cell, group = rcells$group))
alleles <- fread(file.path(prep, 'alleles.tsv.gz'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
res <- run_numbat(count_mat, reference, alleles, genome = genome, t = 1e-5, plot = FALSE,
                  ncores = ncores, out_dir = out)
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
writeLines('Numbat comparator run completed', file.path(out, 'COMPLETE'))
