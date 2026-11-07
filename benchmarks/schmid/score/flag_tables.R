# Add a gate_miss column to every scorer table of one run, written as <name>_gateflag.tsv next to the original
# (originals untouched). gate_miss = TRUE/FALSE for comparator outputs (from gate_methods.tsv), NA for HaploDose and
# for the truth row (not gated). Display names (corr tables) are mapped back with their method_names vector.
# The bin-level merged table (outputs_allmethods_combined*) has no method rows and is not flagged. Tool-level tables
# (filtering_results*) flag a tool if any of its outputs missed. Pairwise corr tables flag method2 and method1 separately.
# Usage: Rscript flag_tables.R GATE_METHODS_TSV RUN_ID EVAL_DIR
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); gm <- fread(a[1], sep = "\t"); rid <- a[2]; e <- a[3]
gm <- gm[run_id == rid]
disp <- c("InferCNV (CNV)" = "infercnv_cnv", "InferCNV (Expr)" = "infercnv_expr", "CaSpER" = "casper", "copyKat" = "copykat",
          "SCEVAN (CNV)" = "scevan_cnv", "SCEVAN (Expr)" = "scevan", "Numbat (Expr)" = "numbat", "Numbat (CNV)" = "numbat_cnv",
          "CONICSmat" = "CONICSmat", "HaploDose" = "haplodose")
# Tool-level rows (filtering_results*): a tool is flagged if any of its outputs missed
tool <- list(inferCNV = c("infercnv_cnv", "infercnv_expr"), SCEVAN = c("scevan_cnv", "scevan"), Numbat = c("numbat", "numbat_cnv"),
             CaSpER = "casper", copyKat = "copykat", CONICSmat = "CONICSmat")
not_gated <- c("haplodose", "HaploDose", "wgs_mean", "WES", "WGS", "CNVkit")
flag_of <- function(nm) {
  key <- ifelse(nm %in% names(disp), disp[nm], nm)
  miss <- gm$gate_miss[match(key, gm$method)]; gap <- gm$max_abs_gap[match(key, gm$method)]
  for (i in which(nm %in% names(tool))) {
    k <- tool[[nm[i]]]; if (!all(k %in% gm$method)) next
    miss[i] <- any(gm$gate_miss[match(k, gm$method)]); gap[i] <- max(gm$max_abs_gap[match(k, gm$method)])
  }
  unk <- setdiff(nm[is.na(miss)], not_gated)
  if (length(unk)) stop("unmapped method names: ", paste(unk, collapse = ", "))
  list(miss = miss, gap = gap)
}
fs <- list.files(e, pattern = "\\.(tsv|txt)$", full.names = TRUE)
fs <- fs[!grepl("_gateflag\\.tsv$|outputs_allmethods_combined|sessionInfo", basename(fs))]
n <- 0
for (f in fs) {
  x <- fread(f, sep = "\t")
  if ("method2" %in% names(x)) {
    # pairwise correlation table: gate_miss refers to method2; method1 gets its own columns
    f2 <- flag_of(x$method2); f1 <- flag_of(x$method1)
    x[, gate_miss := f2$miss]; x[, gate_max_abs_gap := f2$gap]
    x[, gate_miss_method1 := f1$miss]; x[, gate_max_abs_gap_method1 := f1$gap]
  } else if ("method" %in% names(x)) {
    fl <- flag_of(x$method); x[, gate_miss := fl$miss]; x[, gate_max_abs_gap := fl$gap]
  } else { cat("no method column, not flagged:", basename(f), "\n"); next }
  fwrite(x, sub("\\.(tsv|txt)$", "_gateflag.tsv", f), sep = "\t"); n <- n + 1
  cat("flagged", basename(f), "\n")
}
cat("FLAGGED", n, "tables\n")
