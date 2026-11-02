#!/usr/bin/env Rscript
# Kinker 2020 Numbat inputs, stage B: per line of one pool, the steps of numbat's pileup_and_phase.R after the
# pileup (genotype -> Eagle2 per chromosome -> preprocess_allele), on that line's cells only.
# Each channel's stage-A pileup (pileup_channel.sbatch) covers every line in the pool; a line is one person, so its
# cells are taken out of each channel pileup and the per-SNP totals of cellSNP.base.vcf (AD, DP, OTH; used by
# numbat:::genotype and by preprocess_allele's DP_all > 1 & OTH_all == 0 filter) are recomputed from those cells,
# with SNPs below the pileup's --minCOUNT 2 dropped. Cell names are written as Kinker cell ids.
# Output: OUT/<line>/<line>_allele_counts.tsv.gz (channels concatenated), OUT/<line>/{pileup,phasing}/.
# Lines with an existing allele table are skipped (resume).
# args: pileup_root out_dir gmap eagle paneldir ncores channel_title [channel_title ...]
# env PHASE_CHROMS (default 1:22): chromosomes to phase; a subset only for smoke tests on partial pileups.
set.seed(20261102)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 7)
pu_root <- args[1]; outroot <- args[2]; gmap <- args[3]; eagle <- args[4]; paneldir <- args[5]
ncores <- as.integer(args[6]); channels <- args[7:length(args)]
stopifnot(!startsWith(normalizePath(outroot, mustWork = FALSE), '/gpfs/commons/home/jameslee/CRC'))
suppressPackageStartupMessages({
    library(glue); library(stringr); library(data.table); library(dplyr); library(vcfR); library(Matrix)
    library(numbat)
})
chroms <- eval(parse(text = Sys.getenv('PHASE_CHROMS', '1:22')))
genome <- ifelse(str_detect(gmap, 'hg19'), 'hg19', 'hg38')
gtf <- if (genome == 'hg19') gtf_hg19 else gtf_hg38
genetic_map <- fread(gmap) %>% setNames(c('CHROM', 'POS', 'rate', 'cM')) %>% group_by(CHROM) %>%
    mutate(start = POS, end = c(POS[2:length(POS)], POS[length(POS)])) %>% ungroup()

load_channel <- function(ch) {
    d <- file.path(pu_root, ch)
    header <- grep('^##', readLines(file.path(d, 'cellSNP.base.vcf')), value = TRUE)
    fix <- fread(file.path(d, 'cellSNP.base.vcf'), skip = '#CHROM')
    list(name = ch, header = header, fix = fix,
         AD = as(readMM(file.path(d, 'cellSNP.tag.AD.mtx')), 'CsparseMatrix'),
         DP = as(readMM(file.path(d, 'cellSNP.tag.DP.mtx')), 'CsparseMatrix'),
         OTH = as(readMM(file.path(d, 'cellSNP.tag.OTH.mtx')), 'CsparseMatrix'),
         barcodes = fread(file.path(d, 'cellSNP.samples.tsv'), header = FALSE)$V1,
         map = fread(file.path(d, 'cell_map.tsv')))
}
chs <- lapply(channels, load_channel)
lines <- sort(unique(unlist(lapply(chs, function(x) x$map$line))))
message(glue('{length(lines)} lines in channels {paste(channels, collapse = ",")}'))

# Write one channel's pileup restricted to a line's cells, in cellsnp-lite's layout.
subset_channel <- function(x, target, d) {
    # 'target', not 'line': inside data.table's [ ], 'line' would be the column itself (every row matched).
    m <- x$map[x$map$line == target]
    j <- match(m$bam_barcode, x$barcodes)
    j <- j[!is.na(j)]
    if (length(j) == 0) return(FALSE)
    AD <- x$AD[, j, drop = FALSE]; DP <- x$DP[, j, drop = FALSE]; OTH <- x$OTH[, j, drop = FALSE]
    ad <- rowSums(AD); dp <- rowSums(DP); oth <- rowSums(OTH)
    keep <- (dp + oth) >= 2
    if (!any(keep)) return(FALSE)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    fix <- x$fix[keep]
    fix$INFO <- paste0('AD=', ad[keep], ';DP=', dp[keep], ';OTH=', oth[keep])
    writeLines(x$header, file.path(d, 'cellSNP.base.vcf'))
    fwrite(fix, file.path(d, 'cellSNP.base.vcf'), sep = '\t', append = TRUE, col.names = TRUE)
    writeMM(AD[keep, , drop = FALSE], file.path(d, 'cellSNP.tag.AD.mtx'))
    writeMM(DP[keep, , drop = FALSE], file.path(d, 'cellSNP.tag.DP.mtx'))
    fwrite(list(m$cell_id[match(x$barcodes[j], m$bam_barcode)]), file.path(d, 'cellSNP.samples.tsv'), col.names = FALSE)
    TRUE
}

for (line in lines) {
    outdir <- file.path(outroot, line)
    final <- file.path(outdir, glue('{line}_allele_counts.tsv.gz'))
    if (file.exists(final)) { message(glue('skip {line}: exists')); next }
    dir.create(file.path(outdir, 'phasing'), recursive = TRUE, showWarnings = FALSE)
    samples <- c()
    for (x in chs) if (subset_channel(x, line, file.path(outdir, 'pileup', x$name))) samples <- c(samples, x$name)
    if (length(samples) == 0) { message(glue('skip {line}: no SNP reads')); next }
    # From here as pileup_and_phase.R (numbat 1.5.2), with label = line.
    vcfs <- lapply(samples, function(s) vcfR::read.vcfR(file.path(outdir, 'pileup', s, 'cellSNP.base.vcf'), verbose = FALSE))
    vcfs <- lapply(vcfs, function(vcf) { vcf@fix[, 1] <- gsub('chr', '', vcf@fix[, 1]); vcf })
    numbat:::genotype(line, samples, vcfs, file.path(outdir, 'phasing'), chr_prefix = TRUE)
    for (chr in chroms) {
        cmd <- paste(eagle, glue('--numThreads {ncores}'),
                     glue('--vcfTarget {outdir}/phasing/{line}_chr{chr}.vcf.gz'),
                     glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'), glue('--geneticMapFile={gmap}'),
                     glue('--outPrefix {outdir}/phasing/{line}_chr{chr}.phased'))
        status <- system(paste(cmd, '>>', file.path(outdir, 'phasing.log'), '2>&1'))
        if (status != 0) stop(glue('Eagle failed: {line} chr{chr}'))
    }
    vcf_phased <- lapply(chroms, function(chr) {
        fread(glue('{outdir}/phasing/{line}_chr{chr}.phased.vcf.gz'), skip = '#CHROM') %>%
            rename(CHROM = `#CHROM`) %>% mutate(CHROM = str_remove(CHROM, 'chr'))
    }) %>% Reduce(rbind, .) %>% mutate(CHROM = factor(CHROM, unique(CHROM)))
    dfs <- lapply(samples, function(s) {
        pu <- file.path(outdir, 'pileup', s)
        vcf_pu <- fread(file.path(pu, 'cellSNP.base.vcf'), skip = '#CHROM') %>% rename(CHROM = `#CHROM`) %>%
            mutate(CHROM = str_remove(CHROM, 'chr'))
        numbat:::preprocess_allele(sample = line, vcf_pu = vcf_pu, vcf_phased = vcf_phased,
                                   AD = readMM(file.path(pu, 'cellSNP.tag.AD.mtx')),
                                   DP = readMM(file.path(pu, 'cellSNP.tag.DP.mtx')),
                                   barcodes = fread(file.path(pu, 'cellSNP.samples.tsv'), header = FALSE)$V1,
                                   gtf = gtf, gmap = genetic_map) %>% filter(GT %in% c('1|0', '0|1'))
    })
    df <- bind_rows(dfs)
    fwrite(df, paste0(final, '.tmp.gz'), sep = '\t')
    file.rename(paste0(final, '.tmp.gz'), final)
    message(glue('{line}: {nrow(df)} cell-SNP rows, {length(unique(df$cell))} cells, {length(unique(df$snp_id))} SNPs'))
}
capture.output(sessionInfo(), file = file.path(outroot, glue('sessionInfo_{paste(channels, collapse = "_")}.txt')))
