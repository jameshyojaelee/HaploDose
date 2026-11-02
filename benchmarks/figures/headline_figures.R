#!/usr/bin/env Rscript
# Headline benchmark panels for Haplodose. Individual 6 pt PDF panels, Helvetica, useDingbats = FALSE.
# Replot of saved tables only: no statistics are recomputed (setting means and per-line differences are
# arithmetic on the saved per-specimen values). Tests, n and corrections go in captions.md, not on the panels.
#
# Inputs (paths relative to the repo root, ~/Haplodose):
#   A  results/kinker_field_20261102_0529/single_line/rho_by_line.tsv   per-line arm Spearman rho vs DepMap
#      results/kinker_field_20261102_0529/single_line/summary.json      medians, paired Wilcoxon, Holm
#   B  results/kinker_field_20261102_0529/implant/per_specimen.tsv      ARI per line and tool
#      results/kinker_field_20261102_0529/implant/by_setting.tsv        mean ARI per setting and tool
#      results/kinker_field_20261102_0529/implant/summary.json          paired Wilcoxon, Holm over tool pairs
#   C  dev/hd1_proposals_20261102_0600/oracle.tsv                        oracle ARI per line, all region genes
#      dev/hd1_proposals_20261102_0600/oracle_even.tsv                   oracle ARI, Haplodose proposal genes only
#   D  results/andor_rerun_20261102_0042/scores/SNU-601/haplodose/scores.tsv          Haplodose (headline)
#      results/andor_gene_20261102_0113/scores/SNU-601/haplodose/scores.tsv           Haplodose, bin column, v4 run
#      results/andor_gene_20261102_0113/scores/SNU-601/haplodose_gene/scores.tsv      Haplodose, gene column
#      results/andor_20261101_1224/scores_combined/SNU-601/{infercnvpy,infercnv_r,prototype}/scores.tsv
#
# Usage: Rscript benchmarks/figures/headline_figures.R <new_output_dir>   (run from the repo root)
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(jsonlite)})
set.seed(20261102)
source(file.path(Sys.getenv("HOME"), "publication_color_themes.R"))
grDevices::pdf.options(useDingbats = FALSE)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
if (dir.exists(out) && length(list.files(out, pattern = "\\.(pdf|md|tsv)$"))) stop("Use a new output directory")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
FS <- 6
FONT <- "Helvetica"
K <- "results/kinker_field_20261102_0529"
O <- "dev/hd1_proposals_20261102_0600"

theme_hd <- function() {
  theme_classic(base_size = FS, base_family = FONT) %+replace% theme(
    text = element_text(size = FS, family = FONT, face = "plain"),
    axis.text = element_text(size = FS, color = "black"), axis.title = element_text(size = FS),
    legend.text = element_text(size = FS), legend.title = element_blank(),
    strip.text = element_text(size = FS, margin = margin(1, 1, 1, 1)), strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
    legend.key.size = unit(0.22, "cm"), legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 4, 3, 3))
}
save_pdf <- function(p, name, w, h, data = p$data) {
  p <- p + labs(title = NULL, subtitle = NULL, caption = NULL)
  ggsave(file.path(out, name), p, width = w, height = h, device = grDevices::cairo_pdf)
  fwrite(data, file.path(out, sub("\\.pdf$", "_source.tsv", name)), sep = "\t")
}

# Tool display names, colours and shapes. Haplodose = Liang magenta; pseudobulk is the simple baseline and
# takes the control gray; comparators take palette3 stops; older Haplodose variants a light pink.
tools <- data.table(
  id    = c("haplodose", "haplodose_bins", "pseudobulk", "infercnvpy", "infercnv_r", "copykat", "scevan", "prototype"),
  label = c("Haplodose", "Haplodose bins", "Pseudobulk", "infercnvpy", "inferCNV (R)", "CopyKAT", "SCEVAN", "Prototype"),
  color = c("#C9265E", pink_gradient[3], "#9E9E9E", palette3[1], palette3[6], palette3[5], palette3[10], pink_gradient[3]),
  shape = c(16, 1, 18, 17, 2, 15, 4, 1))
col_of <- setNames(tools$color, tools$label)
shp_of <- setNames(tools$shape, tools$label)
lab_of <- setNames(tools$label, tools$id)

# ---- A. Kinker single-line profiles vs DepMap DNA ----------------------------------------------------------
# KEY MESSAGE: Haplodose profiles agree with DNA nearly as well as pseudobulk and better than infercnvpy, CopyKAT, SCEVAN.
rho <- fread(file.path(K, "single_line/rho_by_line.tsv"))
a <- melt(rho, id.vars = "specimen", variable.name = "id", value.name = "rho")
a[, tool := lab_of[as.character(id)]]
ord <- a[, .(m = median(rho)), by = tool][order(m), tool]
a[, tool := factor(tool, levels = ord)]
p <- ggplot(a, aes(rho, tool, color = tool)) +
  geom_jitter(height = 0.22, width = 0, size = 0.35, alpha = 0.55, stroke = 0) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.55, linewidth = 0.25, staplewidth = 0) +
  scale_color_manual(values = col_of, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Spearman ρ vs DNA", y = NULL) + theme_hd()
save_pdf(p, "A1_kinker_singleline_rho.pdf", w = 2.4, h = 1.5)

# KEY MESSAGE: per line, Haplodose is ahead of the published callers; pseudobulk stays slightly ahead of it.
d <- merge(a[id != "haplodose"], rho[, .(specimen, hd = haplodose)], by = "specimen")
d[, diff := hd - rho]
ord <- d[, .(m = median(diff)), by = tool][order(-m), tool]
d[, tool := factor(as.character(tool), levels = ord)]
med <- d[, .(m = median(diff)), by = tool]
p <- ggplot(d, aes(diff, tool, color = tool)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2, color = "#9E9E9E") +
  geom_jitter(height = 0.22, width = 0, size = 0.35, alpha = 0.55, stroke = 0) +
  geom_errorbar(data = med, aes(xmin = m, xmax = m, y = tool), inherit.aes = FALSE, orientation = "y", width = 0.6, linewidth = 0.5, color = "black") +
  scale_color_manual(values = col_of, guide = "none") +
  labs(x = "Δρ, Haplodose − tool", y = NULL) + theme_hd()
save_pdf(p, "A2_kinker_singleline_paired_diff.pdf", w = 2.4, h = 1.3,
         data = d[, .(specimen, tool, rho_haplodose = hd, rho_tool = rho, diff)])

# ---- B. Kinker subclone recovery -------------------------------------------------------------------------
# KEY MESSAGE: all tools recover implanted subclones only for large events; Haplodose and infercnvpy lead.
ps <- fread(file.path(K, "implant/per_specimen.tsv"))
ps[, tool := factor(lab_of[tool], levels = c("Haplodose", "infercnvpy", "CopyKAT", "SCEVAN"))]
ps[, size := factor(size, levels = c(50, 150, 500))]
ps[, direction := factor(direction, levels = c("gain", "loss"), labels = c("Gain", "Loss"))]
ps[, fraction := factor(fraction, levels = c(0.15, 0.4), labels = c("15% of cells", "40% of cells"))]
bs <- ps[, .(mean_ari = mean(ari), n = .N), by = .(tool, size, direction, fraction)]
chk <- fread(file.path(K, "implant/by_setting.tsv"))
stopifnot(isTRUE(all.equal(sort(bs$mean_ari), sort(chk$mean_ari))))  # same means as the saved table
dodge <- position_dodge(width = 0.7)
p <- ggplot(ps, aes(size, ari, color = tool, shape = tool, group = tool)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "#9E9E9E") +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.7, seed = 1),
             size = 0.4, alpha = 0.45, stroke = 0.2) +
  geom_line(data = bs, aes(y = mean_ari), position = dodge, linewidth = 0.35) +
  geom_point(data = bs, aes(y = mean_ari), position = dodge, size = 1.1, stroke = 0.4) +
  facet_grid(direction ~ fraction) +
  scale_color_manual(values = col_of) + scale_shape_manual(values = shp_of) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Event size (genes)", y = "ARI") + theme_hd() +
  theme(legend.position = "top", legend.box.spacing = unit(1, "pt"), panel.spacing = unit(4, "pt"))
save_pdf(p, "B_kinker_subclone_ari.pdf", w = 2.8, h = 2.4)

# ---- C. Headroom: oracle that knows the implanted region vs the tools -------------------------------------
# KEY MESSAGE: the signal is in the counts; even restricted to Haplodose's proposal genes the oracle stays far
# above every tool for events of 50-150 genes.
k <- c("specimen", "size", "direction", "fraction")
orc <- merge(fread(file.path(O, "oracle.tsv"))[, c(k, "oracle_region_em", "oracle_known_params"), with = FALSE],
             fread(file.path(O, "oracle_even.tsv"))[, c(k, "oracle_region_em", "oracle_known_params"), with = FALSE],
             by = k, suffixes = c("", "_proposal"))
stopifnot(nrow(orc) == 166, setequal(orc$specimen, ps$specimen))
os <- orc[, .(all = mean(oracle_region_em), proposal = mean(oracle_region_em_proposal),
              known_all = mean(oracle_known_params), known_proposal = mean(oracle_known_params_proposal), n = .N),
          by = .(size, direction, fraction)]
os[, `:=`(size = factor(size, levels = c(50, 150, 500)),
          direction = factor(direction, levels = c("gain", "loss"), labels = c("Gain", "Loss")),
          fraction = factor(fraction, levels = c(0.15, 0.4), labels = c("15% of cells", "40% of cells")))]
lvl <- os[order(-as.integer(size), -as.integer(direction), -as.integer(fraction)),
          sprintf("%s %s %s", size, tolower(direction), sub(" of cells", "", fraction))]
setting <- function(x) factor(sprintf("%s %s %s", x$size, tolower(x$direction), sub(" of cells", "", x$fraction)), levels = lvl)
bs[, row := setting(bs)]
os[, row := setting(os)]
ol <- melt(os, id.vars = "row", measure.vars = c("all", "proposal"), variable.name = "ref", value.name = "ari")
ol[, ref := factor(ref, levels = c("all", "proposal"), labels = c("Oracle, all genes", "Oracle, proposal genes"))]
p <- ggplot(bs, aes(mean_ari, row)) +
  geom_point(data = ol, aes(ari, row, color = ref), shape = 124, size = 2.6, inherit.aes = FALSE, show.legend = TRUE) +
  geom_point(aes(color = tool, shape = tool), size = 1.0, stroke = 0.4) +
  scale_color_manual(values = c(col_of, "Oracle, all genes" = "#252525", "Oracle, proposal genes" = gray_gradient[8]),
                     breaks = c(levels(ps$tool), levels(ol$ref))) +
  scale_shape_manual(values = shp_of, guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = c(shp_of[levels(ps$tool)], 124, 124), size = c(1, 1, 1, 1, 2.6, 2.6)), ncol = 1)) +
  scale_x_continuous(limits = c(-0.02, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Mean ARI", y = NULL) + theme_hd() + theme(legend.position = "right")
save_pdf(p, "C_kinker_subclone_headroom.pdf", w = 3.0, h = 1.9,
         data = merge(dcast(bs, row ~ tool, value.var = "mean_ari"), os[, .(row, oracle_all = all, oracle_proposal = proposal,
                      oracle_known_all = known_all, oracle_known_proposal = known_proposal, n)], by = "row"))

# ---- D. Andor 2020 SNU-601 with full-depth single-cell DNA truth -------------------------------------------
# KEY MESSAGE: on the one line with full-depth DNA clones, Haplodose finds most clone-differing arms with the
# highest precision; AUROC is comparable across tools.
A <- "results/andor_20261101_1224/scores_combined/SNU-601"
runs <- list(
  haplodose  = "results/andor_rerun_20261102_0042/scores/SNU-601/haplodose/scores.tsv",
  hd_v4_bins = "results/andor_gene_20261102_0113/scores/SNU-601/haplodose/scores.tsv",
  hd_v4_gene = "results/andor_gene_20261102_0113/scores/SNU-601/haplodose_gene/scores.tsv",
  infercnvpy = file.path(A, "infercnvpy/scores.tsv"),
  infercnv_r = file.path(A, "infercnv_r/scores.tsv"),
  prototype  = file.path(A, "prototype/scores.tsv"))
sc <- rbindlist(lapply(names(runs), function(n) fread(runs[[n]])[, .(run = n, source = runs[[n]], diff_precision,
                diff_recall, diff_auroc, diff_called_true, diff_called_arms, dna_subclonal_arms, diff_chance_precision)]))
m <- melt(sc, id.vars = c("run", "source"), measure.vars = c("diff_precision", "diff_recall", "diff_auroc"),
          variable.name = "metric", value.name = "value")
m[, metric := factor(metric, levels = c("diff_precision", "diff_recall", "diff_auroc"), labels = c("Precision", "Recall", "AUROC"))]
rng <- m[run %like% "^h", .(lo = min(value), hi = max(value)), by = metric][, tool := "Haplodose"]
pt <- m[run %in% c("haplodose", "infercnvpy", "infercnv_r", "prototype")]
pt[, tool := factor(lab_of[run], levels = rev(c("Haplodose", "infercnvpy", "inferCNV (R)", "Prototype")))]
rng[, tool := factor(tool, levels = levels(pt$tool))]
chance <- data.table(metric = factor(c("Precision", "AUROC"), levels = levels(m$metric)),
                     x = c(unique(sc$diff_chance_precision), 0.5))
p <- ggplot(pt, aes(value, tool, color = tool, shape = tool)) +
  geom_vline(data = chance, aes(xintercept = x), linewidth = 0.3, linetype = 2, color = "#9E9E9E") +
  geom_errorbar(data = rng, aes(xmin = lo, xmax = hi, y = tool), inherit.aes = FALSE, orientation = "y", width = 0, linewidth = 0.6, color = "#C9265E") +
  geom_point(size = 1.2, stroke = 0.5) + facet_wrap(~metric, nrow = 1) +
  scale_color_manual(values = col_of, guide = "none") + scale_shape_manual(values = shp_of, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = "Clone-differing arms", y = NULL) + theme_hd() + theme(panel.spacing = unit(6, "pt"))
save_pdf(p, "D_andor_snu601_arms.pdf", w = 3.0, h = 1.1, data = sc)

# ---- captions.md: n, unit, test and correction, from the saved summaries ------------------------------------
sj <- fromJSON(file.path(K, "single_line/summary.json"))
pt1 <- as.data.table(sj$paired_tests)
fmt_p <- function(x) formatC(x, format = "g", digits = 2)
lineA <- pt1[, sprintf("%s %d/%d lines higher for Haplodose, Holm p %s", lab_of[b], a_higher, n, fmt_p(p_holm))]
ij <- fromJSON(file.path(K, "implant/summary.json"))
pt2 <- as.data.table(ij$paired_tests)
lineB <- pt2[, sprintf("%s vs %s mean ARI %.3f vs %.3f, Holm p %s", lab_of[a], lab_of[b], mean_a, mean_b, fmt_p(p_holm))]
hd <- sc[run == "haplodose"]
cap_txt <- c(
  "# Headline figure captions", "",
  sprintf("Produced by `benchmarks/figures/headline_figures.R` (seed 20261102). All text 6 pt Helvetica, PDF via cairo_pdf, useDingbats = FALSE. Each panel has a `_source.tsv` with the plotted values."), "",
  "## A1_kinker_singleline_rho.pdf",
  sprintf("Kinker 2020 cell lines profiled one at a time (n = %d lines; unit = cell line). Each point is one line's Spearman rho between the tool's arm-level copy-number profile and DepMap DNA arm copy number; box = median and IQR. Haplodose = gene-level log ratio; Haplodose bins = summed-count bins; pseudobulk = simple baseline. Medians: %s. Source: `%s/single_line/rho_by_line.tsv`, `summary.json`.",
          sj$n_lines, paste(sprintf("%s %.3f", lab_of[names(sj$median_rho)], unlist(sj$median_rho)), collapse = ", "), K), "",
  "## A2_kinker_singleline_paired_diff.pdf",
  sprintf("Per-line difference in rho, Haplodose minus each tool (points right of 0: Haplodose closer to DNA); black tick = median. Paired two-sided Wilcoxon signed-rank per pair, Holm over the 5 pairs with Haplodose, unit = cell line (n = %d): %s. Source: as A1.", sj$n_lines, paste(lineA, collapse = "; ")), "",
  "## B_kinker_subclone_ari.pdf",
  sprintf("Implanted subclonal copy-number events (50, 150 or 500 genes; gain or loss; in 15%% or 40%% of a line's cells), one setting per line (n = %d lines, %d-%d per setting). Small points = adjusted Rand index between a tool's clones and the implanted cells, one per line; large points and lines = setting mean. Overall paired two-sided Wilcoxon, Holm over all 6 tool pairs, unit = cell line: %s. Per-setting differences were not tested. Source: `%s/implant/per_specimen.tsv`, `by_setting.tsv`, `summary.json`.",
          ij$n_specimens, min(bs$n), max(bs$n), paste(lineB, collapse = "; "), K), "",
  "## C_kinker_subclone_headroom.pdf",
  sprintf("Mean ARI per setting for each tool (points, as B) against an oracle that is told the implanted region (black tick: 2-component binomial mixture fitted by EM to each cell's region UMI count, all genes in the region; gray tick: same oracle restricted to the genes Haplodose uses for split proposals). The oracle is a reference bound, not a tool; it ignores overdispersion, so it is optimistic. The oracle with the true parameters known (not drawn) averages %.2f-%.2f per setting with all genes and %.2f-%.2f with proposal genes. Descriptive; no test. n = %d lines, as B. Source: `%s/oracle.tsv`, `%s/oracle_even.tsv` (columns oracle_region_em, oracle_known_params), joined by line to `%s/implant/per_specimen.tsv`.",
          min(os$known_all), max(os$known_all), min(os$known_proposal), max(os$known_proposal), nrow(orc), O, O, K), "",
  "## D_andor_snu601_arms.pdf",
  sprintf("Andor 2020 SNU-601 gastric line, scRNA vs full-depth single-cell DNA (%d of 39 arms differ between some pair of the 15 DNA clones). Precision and recall of arms called as differing between RNA clones, and AUROC of the per-arm clone-difference score. Haplodose: %d/%d called arms true (precision %.2f), recall %d/%d = %.2f, AUROC %.2f; magenta bar = range over three Haplodose scorings of the same code family (the 2026-11-02 rerun, which docs/STATUS.md places on a v6 node, and the gene-column run, which it places on v4, scored on the bin and gene columns), AUROC %.2f-%.2f. Dashed lines: chance precision %.2f, AUROC 0.5. One specimen and 39 arms: no test, no ranking claim. Sources: %s.",
          hd$dna_subclonal_arms, hd$diff_called_true, hd$diff_called_arms, hd$diff_precision, hd$diff_called_true, hd$dna_subclonal_arms,
          hd$diff_recall, hd$diff_auroc, rng[metric == "AUROC", lo], rng[metric == "AUROC", hi], hd$diff_chance_precision,
          paste(sprintf("`%s`", unlist(runs)), collapse = ", ")))
writeLines(cap_txt, file.path(out, "captions.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
