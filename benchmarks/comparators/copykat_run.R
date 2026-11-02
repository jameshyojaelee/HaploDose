#!/usr/bin/env Rscript
# CopyKAT (navinlabcode/copykat @01924c0, env cnv_field), called by run_r_field.py --tool copykat.
# Defaults from the CopyKAT README except norm.cell.names = fit reference cells and output.seg/plot off.
# Writes the CNA matrix (CopyKAT 220 kb genomic windows x cells, smoothed log ratios) as cell-major float32
# with rows.tsv (chromosome, start) and cells.txt, plus calls.tsv (cell_id, malignant, subclone=NA).
# args: prep_dir (as run_infercnv_r.py), out_dir, ncores
set.seed(20261101)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
prep <- args[1]; out <- normalizePath(args[2], mustWork = FALSE); ncores <- as.integer(args[3])
stopifnot(!startsWith(out, '/gpfs/commons/home/jameslee/CRC'))
suppressPackageStartupMessages({library(Matrix); library(copykat)})
dir.create(out, recursive = TRUE, showWarnings = FALSE)
counts <- as(readMM(gzfile(file.path(prep, 'counts.mtx.gz'))), 'CsparseMatrix')
genes <- read.delim(file.path(prep, 'genes.tsv'))$gene
cells <- read.delim(file.path(prep, 'cells.tsv'))$cell_id
rownames(counts) <- genes; colnames(counts) <- cells
ann <- read.delim(file.path(prep, 'annotations.tsv'), header = FALSE)
ref <- ann$V1[grepl('^ref_', ann$V2)]
setwd(out)
t0 <- proc.time()[['elapsed']]
res <- copykat(rawmat = as.matrix(counts), id.type = 'S', cell.line = 'no', ngene.chr = 5, win.size = 25,
               norm.cell.names = ref, KS.cut = 0.1, sam.name = 'run', distance = 'euclidean',
               output.seg = 'FALSE', plot.genes = 'FALSE', genome = 'hg20', n.cores = ncores)
t1 <- proc.time()[['elapsed']]
cna <- res$CNAmat
mat <- as.matrix(cna[, -(1:3)])
key <- setNames(cells, make.names(cells))          # data.frame columns may be make.names()-mangled
cn <- colnames(mat); cn <- ifelse(cn %in% cells, cn, key[cn])
stopifnot(!anyNA(cn))
chrom <- as.character(cna$chrom); chrom[chrom == '23'] <- 'X'
write.table(data.frame(chromosome = paste0('chr', chrom), start = cna$chrompos), file.path(out, 'rows.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
writeLines(cn, file.path(out, 'cells.txt'))
con <- file(file.path(out, 'mat.f32'), 'wb'); writeBin(as.vector(mat), con, size = 4); close(con)  # column-major = cell-major
pred <- res$prediction
pc <- ifelse(pred$cell.names %in% cells, pred$cell.names, key[pred$cell.names])
write.table(data.frame(cell_id = pc, malignant = ifelse(pred$copykat.pred == 'aneuploid', 1,
                                                  ifelse(pred$copykat.pred == 'diploid', 0, NA)), subclone = NA),
            file.path(out, 'calls.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(method = 'copykat', n_cells = ncol(mat), n_rows = nrow(mat), inference_seconds = t1 - t0),
            file.path(out, 'timing.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
