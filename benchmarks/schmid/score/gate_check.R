# Claim 5 gate handling for the scoring job, per STATUS claim 5 (d)/(g): the gate needs at least one primary run whose
# comparator outputs are all within 0.05 of the published Pearson, AUC gains and AUC losses. A comparator output that
# misses on any row is still scored, flagged gate_miss, and its gap reported. A run without a gate table is not scored
# and is reported as missing. No overrides.
# Usage: Rscript gate_check.R INPUTS_TSV REPO_DIR OUT_DIR
# Writes OUT/gate_check.tsv (one row per run x comparator output x metric, plus a row per missing run),
#        OUT/gate_methods.tsv (run_id, method, gate_miss, max_abs_gap) and OUT/runs_to_score.txt.
# Exit 3 if no primary run passes every row.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); map <- fread(a[1], sep = "\t"); H <- a[2]; out <- a[3]
stopifnot(!anyDuplicated(map$run_id))
rows <- list(); meth <- list(); runs <- character(0); primary_pass <- character(0)
for (i in seq_len(nrow(map))) {
  r <- map[i]; f <- file.path(H, r$gate_dir, "gate_comparison.tsv")
  if (!file.exists(f) || file.size(f) == 0) {
    rows[[length(rows) + 1]] <- data.table(run_id = r$run_id, dataset = r$dataset, role = r$role, status = "gate_table_missing",
      method = NA, metric = NA, ours = NA, published = NA, difference = NA, gate_miss = NA, gate_table_md5 = NA)
    next
  }
  g <- fread(f, sep = "\t"); md5 <- unname(tools::md5sum(f))
  stopifnot(all(c("method", "metric", "ours", "published", "difference") %in% names(g)), nrow(g) > 0)
  g[, row_miss := is.na(difference) | abs(difference) > 0.05]
  rows[[length(rows) + 1]] <- data.table(run_id = r$run_id, dataset = r$dataset, role = r$role, status = "gate_table_present",
    method = g$method, metric = g$metric, ours = g$ours, published = g$published, difference = g$difference,
    gate_miss = g$row_miss, gate_table_md5 = md5)
  # A method with no output (all differences NA, e.g. Numbat on MCF7) is a gate miss with max_abs_gap NA.
  m <- g[, .(gate_miss = any(row_miss),
             max_abs_gap = if (all(is.na(difference))) NA_real_ else max(abs(difference), na.rm = TRUE)), by = method]
  meth[[length(meth) + 1]] <- cbind(run_id = r$run_id, m)
  runs <- c(runs, r$run_id)
  if (r$role == "primary" && !any(g$row_miss)) primary_pass <- c(primary_pass, r$run_id)
}
fwrite(rbindlist(rows), file.path(out, "gate_check.tsv"), sep = "\t")
fwrite(if (length(meth)) rbindlist(meth) else data.table(run_id = character(), method = character(), gate_miss = logical(),
       max_abs_gap = numeric()), file.path(out, "gate_methods.tsv"), sep = "\t")
writeLines(runs, file.path(out, "runs_to_score.txt"))
writeLines(primary_pass, file.path(out, "primary_runs_passing.txt"))
cat("runs with gate table:", runs, "\nprimary runs passing every row:", primary_pass, "\n")
if (!length(primary_pass)) { cat("NO_PRIMARY_GATE_PASS\n"); quit(status = 3) }
cat("GATE_CHECK_OK\n")
