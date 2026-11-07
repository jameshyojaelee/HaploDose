# Scorer fixture (no HaploDose output used): a fake "HaploDose" output built from the gate's copyKat output, so the
# patched scorers must give HaploDose exactly copyKat's metrics. One clone fx_c1 = copyKat pseudobulk log ratio over
# non-reference cells; every non-reference copyKat cell is assigned to it.
# Usage: Rscript make_fixture.R GATE_PIPELINE_DIR DS OUT_HD_DIR
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
a <- commandArgs(TRUE); gp <- a[1]; ds <- a[2]; out <- a[3]
source(file.path(gp, "scripts/lib.R"))
inp <- file.path(gp, "data", paste0("input_", ds))
ck <- read_copykat(file.path(gp, "results", paste0("output_", ds), "copykat", paste0(ds, "_copykat_CNA_raw_results_gene_by_cell.txt")),
                   file.path(inp, "sample_annotation.txt"), file.path(inp, "ref_groups.txt"))
g <- ck[[1]]
dir.create(file.path(out, "haplodose"), recursive = TRUE)
fwrite(data.frame(gene = paste0("g", seq_along(g)), chromosome = as.character(seqnames(g)), start = start(g), end = end(g),
                  fx_c1 = g$copykat - 2), file.path(out, "haplodose", paste0(ds, "_haplodose_gene_by_clone.tsv")), sep = "\t")
cr <- fread(file.path(gp, "results", paste0("output_", ds), "copykat", paste0(ds, "_copykat_CNA_raw_results_gene_by_cell.txt")), nrows = 1)
an <- fread(file.path(inp, "sample_annotation.txt"), header = FALSE); rg <- fread(file.path(inp, "ref_groups.txt"))
cells <- intersect(colnames(cr)[8:ncol(cr)], an$V1[!(an$V2 %in% rg$ref_groups)])
stopifnot(length(cells) == ck[[2]])
fwrite(data.frame(cell = cells, clone = "fx_c1"), file.path(out, "haplodose", paste0(ds, "_haplodose_cell_clone.tsv")), sep = "\t")
cat("FIXTURE", ds, length(g), "genes", length(cells), "cells\n")
