#!/usr/bin/env Rscript
# SCEVAN (AntonioDeFalco/SCEVAN @5a49b88, v1.0.3, env cnv_field), called by run_r_field.py --tool scevan.
# pipelineCNA with norm_cell = fit reference cells and SUBCLONES = TRUE; other settings default.
# Writes SCEVAN's relative CNA matrix (genes x tumour cells, from output/run_CNAmtx.RData) as cell-major
# float32 with rows.tsv (chromosome, start) and cells.txt, plus calls.tsv (cell_id, malignant, subclone).
# args: prep_dir (as run_infercnv_r.py), out_dir, ncores
set.seed(20261101)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
prep <- args[1]; out <- normalizePath(args[2], mustWork = FALSE); ncores <- as.integer(args[3])
stopifnot(!startsWith(out, '/gpfs/commons/home/jameslee/CRC'))
suppressPackageStartupMessages({library(Matrix); library(SCEVAN)})
# SCEVAN's plotConsensusCNA indexes dfL[[grep(paste0('subclone', x), names(dfL))]], which matches subclone1
# and subclone10.. when there are >= 10 subclones and stops ("no such index at level 2"). It only draws
# consensus.png, so it is replaced with a no-op; no output used here depends on it.
assignInNamespace('plotConsensusCNA', function(...) invisible(NULL), ns = 'SCEVAN')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
counts <- as(readMM(gzfile(file.path(prep, 'counts.mtx.gz'))), 'CsparseMatrix')
genes <- read.delim(file.path(prep, 'genes.tsv'))$gene
cells <- read.delim(file.path(prep, 'cells.tsv'))$cell_id
rownames(counts) <- genes; colnames(counts) <- cells
ann <- read.delim(file.path(prep, 'annotations.tsv'), header = FALSE)
ref <- ann$V1[grepl('^ref_', ann$V2)]
setwd(out)
t0 <- proc.time()[['elapsed']]
res <- pipelineCNA(as.matrix(counts), sample = 'run', par_cores = ncores, norm_cell = ref, SUBCLONES = TRUE,
                   plotTree = FALSE, organism = 'human')
t1 <- proc.time()[['elapsed']]
e <- new.env(); load(file.path(out, 'output', 'run_CNAmtx.RData'), envir = e)
print(ls(e))
mat <- as.matrix(get(ls(e)[1], envir = e))
a <- new.env(); load(file.path(out, 'output', 'run_count_mtx_annot.RData'), envir = a)
print(ls(a))
annot <- get(ls(a)[1], envir = a)
print(head(annot)); print(dim(mat)); print(dim(annot))
stopifnot(nrow(annot) == nrow(mat))
chrom <- as.character(annot$seqnames); chrom <- sub('^chr', '', chrom); chrom[chrom == '23'] <- 'X'
write.table(data.frame(chromosome = paste0('chr', chrom), start = annot$start), file.path(out, 'rows.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
writeLines(colnames(mat), file.path(out, 'cells.txt'))
con <- file(file.path(out, 'mat.f32'), 'wb'); writeBin(as.vector(mat), con, size = 4); close(con)
sub <- if ('subclone' %in% names(res)) res$subclone else NA
write.table(data.frame(cell_id = rownames(res), malignant = ifelse(res$class == 'tumor', 1, ifelse(res$class == 'normal', 0, NA)),
                       subclone = sub), file.path(out, 'calls.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(data.frame(method = 'scevan', n_cells = ncol(mat), n_rows = nrow(mat), inference_seconds = t1 - t0),
            file.path(out, 'timing.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
