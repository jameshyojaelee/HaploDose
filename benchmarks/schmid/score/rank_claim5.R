# Claim 5 ranking (STATUS item 5 and amendment a-c, 2026-11-02; missing-output rule 2026-11-04 ~16:45).
# Per dataset and metric (Pearson vs DNA truth, AUC gains, AUC losses; higher is better) the 10 outputs are ranked
# 1 = best. Ties share the average rank. An output missing for a dataset and metric ranks last there; several missing
# outputs share the last ranks (average). SOTA = HaploDose alone has the best (lowest) mean rank over
# 3 metrics x the primary datasets. Descriptive: the paper's normalised overall score (mean of Pearson floored at 0,
# (AUC - 0.5) x 2 floored at 0 for gains and losses, and max F1), averaged over the same datasets.
# Inputs: a TSV with columns set, run_id, eval_dir (the evaluation_haplodose dir of one scored run). Rows are grouped by
# set (e.g. primary, sens_mm568, permethod...). A set's rows must be distinct datasets.
# Usage: Rscript rank_claim5.R RUNS_TSV OUT_DIR [MODE]   MODE = combined (default) or permethod
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); runs <- fread(a[1], sep = "\t"); out <- a[2]; mode <- if (length(a) >= 3) a[3] else "combined"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
keys <- c("infercnv_cnv", "infercnv_expr", "casper", "copykat", "scevan_cnv", "scevan", "numbat", "numbat_cnv",
          "CONICSmat", "haplodose")
disp <- c("InferCNV (CNV)" = "infercnv_cnv", "InferCNV (Expr)" = "infercnv_expr", "CaSpER" = "casper",
          "copyKat" = "copykat", "SCEVAN (CNV)" = "scevan_cnv", "SCEVAN (Expr)" = "scevan", "Numbat (Expr)" = "numbat",
          "Numbat (CNV)" = "numbat_cnv", "CONICSmat" = "CONICSmat", "HaploDose" = "haplodose")
one <- function(d) {
  if (mode == "combined") {
    fc <- list.files(d, "^evaluation_cnv_prediction_corr.*_gateflag\\.tsv$", full.names = TRUE)
    fa <- list.files(d, "^evaluation_cnv_prediction_auc.*_gateflag\\.tsv$", full.names = TRUE)
    ff <- list.files(d, "^evaluation_cnv_prediction_f1.*_gateflag\\.tsv$", full.names = TRUE)
    stopifnot(length(fc) == 1, length(fa) == 1, length(ff) == 1)
    cr <- fread(fc); cr <- cr[method1 %in% c("WES", "WGS", "CNVkit") & method2 %in% names(disp)]
    p <- data.table(method = unname(disp[cr$method2]), pearson = cr$pearson)
    au <- fread(fa)[, .(method, auc_gains, auc_losses)]
    f1 <- fread(ff)[, .(method, max_f1)]
    x <- merge(merge(p, au, by = "method", all = TRUE), f1, by = "method", all = TRUE)
  } else {
    fp <- list.files(d, "^evaluation_cnv_prediction_permethod.*_gateflag\\.tsv$", full.names = TRUE)
    stopifnot(length(fp) == 1)
    x <- fread(fp)[, .(method, pearson, auc_gains, auc_losses, max_f1)]
  }
  stopifnot(!anyDuplicated(x$method))
  unk <- setdiff(x$method, keys); if (length(unk)) stop("unknown method keys: ", paste(unk, collapse = ", "))
  x <- merge(data.table(method = keys), x, by = "method", all.x = TRUE)
  x
}
rk <- function(v) {  # higher is better; NA ranks last, ties averaged (rank() with na.last = "keep" then fill)
  r <- rank(-v, ties.method = "average", na.last = "keep")
  nmiss <- sum(is.na(v)); n <- length(v)
  if (nmiss) r[is.na(v)] <- mean((n - nmiss + 1):n)
  r
}
long <- rbindlist(lapply(seq_len(nrow(runs)), function(i) cbind(set = runs$set[i], run_id = runs$run_id[i], one(runs$eval_dir[i]))))
for (s in unique(long$set)) stopifnot(!anyDuplicated(runs[set == s]$run_id))
long[, `:=`(rank_pearson = rk(pearson), rank_auc_gains = rk(auc_gains), rank_auc_losses = rk(auc_losses)), by = .(set, run_id)]
long[, norm_score := rowMeans(cbind(pmax(pearson, 0), pmax((auc_gains - 0.5) * 2, 0), pmax((auc_losses - 0.5) * 2, 0), max_f1))]
fwrite(long, file.path(out, paste0("per_dataset_", mode, ".tsv")), sep = "\t")
summ <- long[, .(n_datasets = .N, mean_rank = mean(c(rank_pearson, rank_auc_gains, rank_auc_losses)),
                 n_missing_cells = sum(is.na(pearson)) + sum(is.na(auc_gains)) + sum(is.na(auc_losses)),
                 mean_norm_score = mean(norm_score)), by = .(set, method)]
summ[, order_by_mean_rank := rank(mean_rank, ties.method = "min"), by = set]
setorder(summ, set, mean_rank)
fwrite(summ, file.path(out, paste0("summary_", mode, ".tsv")), sep = "\t")
verdict <- summ[, {
  hd <- mean_rank[method == "haplodose"]; best <- min(mean_rank[method != "haplodose"])
  .(haplodose_mean_rank = hd, best_other_mean_rank = best, best_other = paste(method[method != "haplodose" & mean_rank == best], collapse = ","),
    sota = hd < best)
}, by = set]
fwrite(verdict, file.path(out, paste0("verdict_", mode, ".tsv")), sep = "\t")
print(verdict)
cat("RANK_DONE", mode, "\n")
