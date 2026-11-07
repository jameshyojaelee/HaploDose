# Claim 5 deviation (docs/STATUS.md, claim 5, 2026-11-04): build the two gene-annotation tables that their pipeline
# fetches live with biomaRt (CaSpER run_casper_gene_annot, CONICSmat CONICSmat_import_data), from raw BioMart tables
# downloaded once from the Ensembl 113 archive (oct2024.archive.ensembl.org) by fetch_biomart.sh. biomaRt 2.50.3 in
# schmid_r41 cannot query Ensembl any more (HTTP 404 / BiocFileCache-dbplyr error).
# getBM(filters = "hgnc_symbol", values = genes) is replaced by subsetting the unfiltered raw table to
# hgnc_symbol %in% genes, keeping the server's row order. Everything after getBM is copied verbatim from
# CaSpER::generateAnnotation (src/CaSpER/R/utility_functions.R, commit 8aae116) and CONICSmat::getGenePositions
# (src/CONICS/CONICSmat/R/GetPositions.R, commit 053cf26), and the write.table calls are those of
# casper_gene_annotations.R and CONICSmat_import_data.R.
# Usage: Rscript build_annotations.R RAW_DIR COUNT_MATRIX DATASET OUT_DIR
suppressMessages({library(data.table); library(CaSpER)})
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1]; input_file <- args[2]; dataset <- args[3]; out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_raw <- function(f) {
  x <- readLines(f)
  stopifnot(tail(x, 1) == "[success]")          # BioMart completionStamp: the table is complete
  fread(text = paste(head(x, -1), collapse = "\n"), sep = "\t", header = TRUE, colClasses = "character",
        na.strings = NULL, data.table = FALSE)
}

# ---- CaSpER (their casper_gene_annotations.R: genes = count_matrix$V1 from fread) ----
count_matrix <- fread(input_file, select = 1)
genes <- count_matrix[[1]]
data(hg38_cytoband)
raw_c <- read_raw(file.path(raw_dir, "biomart_casper.tsv"))
colnames(raw_c) <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band")
raw_c$start_position <- as.integer(raw_c$start_position); raw_c$end_position <- as.integer(raw_c$end_position)
# getBM attributes were c(id_type, "hgnc_symbol", ...) with id_type = "hgnc_symbol": the symbol column twice
G_list <- raw_c[raw_c$hgnc_symbol %in% genes, c(1, 1:5)]
colnames(G_list) <- c("hgnc_symbol", "hgnc_symbol.1", "chromosome_name", "start_position", "end_position", "band")
id_type <- "hgnc_symbol"; centromere <- centromere_hg38
# --- verbatim from generateAnnotation, after getBM ---
common <- intersect(genes, G_list[, id_type])
ord <- match(common, G_list[, id_type])
annotation <- G_list[ord, ]
annotation <- annotation[order(annotation$start_position), ]
annotation$cytoband <- paste0(annotation$chromosome_name, substr((annotation$band), 0, 1))
chr.list <- as.character(1:22, "X")
idx <- unlist(unique(as.vector(sapply(chr.list, function(x) as.vector(unlist(which(as.character(annotation$chromosome_name) ==
    x)))))))
annotation <- as.data.frame(annotation[idx, ])
colnames(annotation)[c(1:5, 7)] <- c("Gene", "GeneSymbol",  "Chr", "start", "end", "cytoband")
annotation$isCentromer <- rep("no", nrow(annotation))
centromere_snps <- NULL
for (k in 1:(dim(centromere)[1])) {
    annotation$isCentromer[which(as.character(annotation$Chr) == gsub("chr", "", as.character(centromere$V1[k])) & (as.numeric(as.character(annotation$Position)) >=
        centromere$V2[k] & as.numeric(as.character(annotation$Position)) <= centromere$V3[k]))] <- "yes"
}
annotation$Position <- (as.numeric(annotation$start) + as.numeric(annotation$end))/2
annotation$new_positions <- as.vector(unlist(lapply(lapply(split(annotation$cytoband, annotation$cytoband), length)[unique(annotation$cytoband)],
    function(x) 1:x)))
# --- end verbatim ---
write.table(annotation, file = file.path(out_dir, paste0("gene_annotation_casper_", dataset, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- CONICSmat (their CONICSmat_import_data.R: genes = rownames of read.table(check.names = F)) ----
gene_names <- genes   # same first column; read.table(row.names = 1) does not alter row names
stopifnot(!anyDuplicated(gene_names))
raw_k <- read_raw(file.path(raw_dir, "biomart_conicsmat.tsv"))
colnames(raw_k) <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")
raw_k$start_position <- as.integer(raw_k$start_position); raw_k$end_position <- as.integer(raw_k$end_position)
gene_positions <- raw_k[raw_k$hgnc_symbol %in% gene_names, ]
ignoreAlt <- F
# --- verbatim from getGenePositions, after getBM ---
if (ignoreAlt == T){
gene_positions=gene_positions[order(gene_positions[,3],decreasing=F),]
	}
gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
gene_positions[which(gene_positions[,3]=="X"),3]=23
gene_positions[which(gene_positions[,3]=="Y"),3]=24
gene_positions[which(gene_positions[,3]=="MT"),3]=0
gene_positions[which(nchar(gene_positions[,3])>2),3]=0
gene_positions=gene_positions[order(as.numeric(gene_positions[,3]),decreasing=F),]
# --- end verbatim ---
write.table(gene_positions, file = file.path(out_dir, paste0("gene_annotation_conicsmat_", dataset, ".tsv")),
            sep = "\t", row.names = FALSE, eol = "\r")
cat(dataset, "genes", length(genes), "casper_rows", nrow(annotation), "conicsmat_rows", nrow(gene_positions), "\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
