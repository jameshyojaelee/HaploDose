# ------------------------------------------------------------------------------
# HaploDose readers for the Schmid et al. 2025 benchmark scorer
# (external/benchmark_scrnaseq_cnv_callers/snakemake_pipeline/scripts/lib.R).
# Source after their lib.R. Mirrors read_copykat / read_copykat_individual_cells:
# returns list(per-gene GRanges with a score column where 2 = neutral, number of cells).
#
# Inputs are written by benchmarks/schmid/run_haplodose.py:
#   gene_clone_path  <dataset>_haplodose_gene_by_clone.tsv: gene, chromosome, start, end, one column per
#                    clone (natural log ratio, 0 = neutral, the value of the bin that holds the gene)
#   cell_clone_path  <dataset>_haplodose_cell_clone.tsv: cell, clone (query cells with an assigned clone)
# Each non-reference cell takes its MAP clone's profile, so the pseudobulk is the cell-count-weighted mean of
# the clone profiles over the non-reference cells HaploDose assigned. As for copyKat (log ratio + 2), the score
# is the log ratio + 2. Pearson and both AUCs do not change under this shift; the truncated AUC threshold at 2
# then sits at log ratio 0.
# ------------------------------------------------------------------------------

read_haplodose<-function(gene_clone_path,cell_clone_path,annot_path,ref_path,clones=NULL){

  hd<-fread(gene_clone_path)
  hd_grange<-makeGRangesFromDataFrame(hd,seqnames.field="chromosome",
                                      start.field="start",end.field="end")
  clone_cols<-setdiff(colnames(hd),c("gene","chromosome","start","end"))

  #Remove reference cells (same rule as read_copykat)
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  cells<-fread(cell_clone_path)
  cells<-cells[cells$cell %in% cancer_cells & cells$clone %in% clone_cols,]
  if(nrow(cells)==0){
    stop("No non-reference cells with an assigned HaploDose clone")
  }

  if(is.null(clones)){
    n_per_clone<-table(factor(cells$clone,levels=clone_cols))
    w<-as.numeric(n_per_clone)/sum(n_per_clone)
    hd_grange$haplodose<-as.vector(as.matrix(hd[,..clone_cols]) %*% w)+2
  } else {
    #Per-clone pseudobulk (columns clone_<name>), as read_copykat does when a clone table is given
    pbulk_clones<-as.matrix(hd[,..clone_cols])+2
    colnames(pbulk_clones)<-paste0("clone_",clone_cols)
    mcols(hd_grange)<-cbind(mcols(hd_grange),pbulk_clones)
  }

  return(list(hd_grange,nrow(cells)))
}

#Per-cell matrix (genes x non-reference assigned cells), each cell holding its MAP clone's profile + 2
read_haplodose_individual_cells<-function(gene_clone_path,cell_clone_path,annot_path,ref_path){

  hd<-fread(gene_clone_path)
  hd_grange<-makeGRangesFromDataFrame(hd,seqnames.field="chromosome",
                                      start.field="start",end.field="end")
  clone_cols<-setdiff(colnames(hd),c("gene","chromosome","start","end"))

  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  cells<-fread(cell_clone_path)
  cells<-cells[cells$cell %in% cancer_cells & cells$clone %in% clone_cols,]

  hd_cells<-as.matrix(hd[,..clone_cols])[,cells$clone,drop=FALSE]+2
  colnames(hd_cells)<-cells$cell

  return(list(hd_grange,hd_cells))
}
