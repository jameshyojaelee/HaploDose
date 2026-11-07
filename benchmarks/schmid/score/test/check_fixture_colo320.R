# Checks for run_fixture_colo320.sbatch. Usage: Rscript check_fixture_colo320.R SCORE_OUT RC
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); A <- a[1]; rc <- a[2]; ds <- "COLO320"
ok <- TRUE; chk <- function(x, msg) { cat(if (isTRUE(x)) "ok  " else "FAIL", msg, "\n"); if (!isTRUE(x)) ok <<- FALSE }
chk(rc == "0", paste("run_score_all exit", rc))
gm <- fread(file.path(A, "gate_methods.tsv"))
chk(isTRUE(gm[run_id == ds & method == "casper"]$gate_miss), "gate_methods: casper gate_miss TRUE")
chk(!isTRUE(gm[run_id == ds & method == "copykat"]$gate_miss), "gate_methods: copykat gate_miss FALSE")
sk <- fread(file.path(A, "skipped_runs.tsv"))
chk(all(sk$run_id == "BCC06post_gateanchor"), "only the anchor is skipped")
for (arm in c("vnext", "sens1424")) {
  e <- file.path(A, arm, ds, "pipeline", "results", paste0("output_", ds), "evaluation_haplodose")
  chk(file.exists(file.path(e, "evaluation_cnv_prediction_auc_cnvkit.txt")), paste(arm, ": combined AUC table (cnvkit) written"))
  chk(file.exists(file.path(e, "evaluation_cnv_prediction_permethod_cnvkit.tsv")), paste(arm, ": per-method table written"))
  fl <- list.files(e, pattern = "_gateflag\\.tsv$", full.names = TRUE)
  chk(length(fl) >= 4, paste(arm, ": flagged tables written:", length(fl)))
  for (f in fl) {
    x <- fread(f); mc <- intersect(c("method2", "method"), names(x))[1]
    if (mc == "method2") x <- x[method1 %in% c("WES", "WGS", "CNVkit")]
    key <- tolower(x[[mc]]); h <- x[key == "haplodose"]; k <- x[key == "copykat"]; cs <- x[key == "casper"]
    if (nrow(h) != 1 || nrow(k) != 1) { chk(FALSE, paste(arm, basename(f), "rows haplodose", nrow(h), "copykat", nrow(k))); next }
    num <- setdiff(names(x)[sapply(x, is.numeric)], c("gate_max_abs_gap", "gate_max_abs_gap_method1"))
    dif <- max(abs(unlist(h[, ..num]) - unlist(k[, ..num])), na.rm = TRUE)
    chk(dif < 1e-9, sprintf("%s %s: max|haplodose-copykat| = %.3g", arm, basename(f), dif))
    chk(is.na(h$gate_miss), paste(arm, basename(f), ": HaploDose gate_miss NA"))
    chk(identical(k$gate_miss, FALSE), paste(arm, basename(f), ": copyKat gate_miss", k$gate_miss))
    if (nrow(cs) == 1) chk(identical(cs$gate_miss, TRUE), paste(arm, basename(f), ": CaSpER gate_miss", cs$gate_miss))
  }
}
cat(if (ok) "FIXTURE_PASS" else "FIXTURE_FAIL", "\n")
