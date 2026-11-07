# Checks for run_fixture.sbatch cases A, B and (if DS2 != none) C. Usage: Rscript check_fixture.R OUT DS EXIT_A EXIT_B EXIT_C DS2
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); d <- a[1]; ds <- a[2]; rcA <- as.integer(a[3]); rcB <- as.integer(a[4]); ok <- TRUE
rcC <- suppressWarnings(as.integer(a[5])); ds2 <- if (length(a) >= 6) a[6] else "none"
chk <- function(cond, msg) { cat(if (cond) "ok  " else "FAIL", msg, "\n"); if (!cond) ok <<- FALSE }
A <- file.path(d, "score_A"); B <- file.path(d, "score_B")
chk(rcA == 0, "case A exits 0")
chk(rcB == 3, "case B (no primary passes) exits 3")
chk(grepl("GATE_CHECK_STOP", paste(readLines(file.path(d, "case_B.out")), collapse = "\n")), "case B prints GATE_CHECK_STOP")
chk(!dir.exists(file.path(B, "vnext")) && !dir.exists(file.path(B, "sens1424")), "case B scored nothing")
gc <- fread(file.path(A, "gate_check.tsv")); gm <- fread(file.path(A, "gate_methods.tsv"))
sk <- fread(file.path(A, "skipped_runs.tsv"))
chk(setequal(sk$run_id, paste0(ds, c("_nogate", "_badinput"))) && nrow(sk) == 2, "case A: exactly nogate and badinput skipped")
chk(sk[run_id == paste0(ds, "_nogate")]$reason == "gate_table_missing", "case A: nogate reason gate_table_missing")
chk(grepl("^haplodose_input_not_identical_to_gate:", sk[run_id == paste0(ds, "_badinput")]$reason), "case A: badinput reason input not identical")
chk(nrow(gc[run_id == paste0(ds, "_nogate") & status == "gate_table_missing"]) == 1, "case A: gate_check.tsv reports the missing table")
chk(isTRUE(gm[run_id == paste0(ds, "_miss") & method == "copykat"]$gate_miss) &&
    abs(gm[run_id == paste0(ds, "_miss") & method == "copykat"]$max_abs_gap - 0.07) < 1e-12, "case A: copykat gate_miss TRUE with gap 0.07 (miss run)")
chk(!any(gm[run_id == paste0(ds, "_miss") & method != "copykat"]$gate_miss), "case A: other methods not flagged (miss run)")
chk(!any(gm[run_id == ds]$gate_miss), "case A: real run has no gate_miss")
nn <- gm[run_id == paste0(ds, "_nonumbat")]
chk(nrow(nn[method %in% c("numbat", "numbat_cnv") & gate_miss == TRUE & is.na(max_abs_gap)]) == 2,
    "case A: nonumbat run has numbat and numbat_cnv gate_miss TRUE with max_abs_gap NA")
chk(!any(nn[!(method %in% c("numbat", "numbat_cnv"))]$gate_miss) && all(is.finite(nn[!(method %in% c("numbat", "numbat_cnv"))]$max_abs_gap)),
    "case A: nonumbat run, other methods not flagged and finite gaps")
chk(nrow(gc[run_id == paste0(ds, "_nonumbat") & method %in% c("numbat", "numbat_cnv") & is.na(ours) & gate_miss == TRUE]) == 6,
    "case A: gate_check.tsv has the 6 Numbat rows with ours NA and gate_miss TRUE")
for (arm in c("vnext", "sens1424")) {
  chk(!dir.exists(file.path(A, arm, paste0(ds, "_nogate"))), paste("case A", arm, ": nogate run not scored"))
  chk(!dir.exists(file.path(A, arm, paste0(ds, "_badinput"))), paste("case A", arm, ": badinput run not scored"))
  for (rid in c(ds, paste0(ds, "_miss"), paste0(ds, "_nonumbat"))) {
    e <- file.path(A, arm, rid, "pipeline", "results", paste0("output_", ds), "evaluation_haplodose")
    fl <- list.files(e, pattern = "_gateflag\\.tsv$", full.names = TRUE)
    chk(length(fl) >= 4, paste("case A", arm, rid, ": flagged tables written:", length(fl)))
    for (f in fl) {
      x <- fread(f); mc <- intersect(c("method2", "method"), names(x))[1]
      if (mc == "method2") x <- x[method1 %in% c("WES", "WGS", "CNVkit")]   # truth rows of the pairwise table
      key <- tolower(x[[mc]])
      h <- x[key == "haplodose"]; k <- x[key == "copykat"]
      if (nrow(h) != 1 || nrow(k) != 1) { chk(FALSE, paste(basename(f), "rows haplodose", nrow(h), "copykat", nrow(k))); next }
      num <- setdiff(names(x)[sapply(x, is.numeric)], c("gate_max_abs_gap", "gate_max_abs_gap_method1"))
      dif <- max(abs(unlist(h[, ..num]) - unlist(k[, ..num])), na.rm = TRUE)
      chk(dif < 1e-9, sprintf("%s %s %s: max|haplodose-copykat| = %.3g", arm, rid, basename(f), dif))
      chk(is.na(h$gate_miss), paste(arm, rid, basename(f), ": HaploDose gate_miss NA"))
      chk(identical(k$gate_miss, rid == paste0(ds, "_miss")), paste(arm, rid, basename(f), ": copyKat gate_miss", k$gate_miss))
      isnb <- grepl("^numbat", key)
      if (rid == paste0(ds, "_nonumbat")) {
        nb <- x[isnb]
        chk(nrow(nb) >= 1 && all(nb$gate_miss %in% TRUE), paste(arm, rid, basename(f), ": Numbat rows gate_miss TRUE:", nrow(nb)))
      }
      oth <- x[!(key %in% c("haplodose", "copykat")) & !(isnb & rid == paste0(ds, "_nonumbat")) & !is.na(gate_miss)]
      chk(!any(oth$gate_miss), paste(arm, rid, basename(f), ": other comparators gate_miss FALSE"))
    }
  }
}
if (ds2 != "none") {
  C <- file.path(d, "score_C")
  chk(identical(rcC, 0L), "case C exits 0")
  gmC <- fread(file.path(C, "gate_methods.tsv"))
  chk(nrow(gmC[run_id == ds2 & method %in% c("numbat", "numbat_cnv") & gate_miss == TRUE & is.na(max_abs_gap)]) == 2,
      "case C: DS2 numbat and numbat_cnv gate_miss TRUE with max_abs_gap NA")
  for (arm in c("vnext", "sens1424")) {
    e <- file.path(C, arm, ds2, "pipeline", "results", paste0("output_", ds2), "evaluation_haplodose")
    chk(file.exists(file.path(e, "evaluation_cnv_prediction_auc.txt")), paste("case C", arm, ": combined AUC table written"))
    chk(file.exists(file.path(e, "evaluation_cnv_prediction_permethod.tsv")), paste("case C", arm, ": per-method table written"))
    fl <- list.files(e, pattern = "_gateflag\\.tsv$", full.names = TRUE)
    chk(length(fl) >= 4, paste("case C", arm, ": flagged tables written:", length(fl)))
    for (f in fl) {
      x <- fread(f); mc <- intersect(c("method2", "method"), names(x))[1]
      allkeys <- tolower(c(x[[mc]], if (mc == "method2") x$method1))
      chk(!any(grepl("numbat", allkeys)), paste("case C", arm, basename(f), ": no Numbat rows"))
      if (mc == "method2") x <- x[method1 %in% c("WES", "WGS", "CNVkit")]
      key <- tolower(x[[mc]]); h <- x[key == "haplodose"]; k <- x[key == "copykat"]
      if (nrow(h) != 1 || nrow(k) != 1) { chk(FALSE, paste("case C", basename(f), "rows haplodose", nrow(h), "copykat", nrow(k))); next }
      num <- setdiff(names(x)[sapply(x, is.numeric)], c("gate_max_abs_gap", "gate_max_abs_gap_method1"))
      dif <- max(abs(unlist(h[, ..num]) - unlist(k[, ..num])), na.rm = TRUE)
      chk(dif < 1e-9, sprintf("case C %s %s: max|haplodose-copykat| = %.3g", arm, basename(f), dif))
      chk(is.na(h$gate_miss), paste("case C", arm, basename(f), ": HaploDose gate_miss NA"))
      chk(identical(k$gate_miss, gmC[run_id == ds2 & method == "copykat"]$gate_miss), paste("case C", arm, basename(f), ": copyKat flag matches gate_methods"))
    }
  }
}
cat(if (ok) "FIXTURE_PASS\n" else "FIXTURE_FAIL\n")
