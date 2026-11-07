# Write count_matrix.txt exactly as the Schmid preprocessing scripts do (write.table, sep tab, quote FALSE,
# row.names TRUE), from the counts.mtx / genes.txt / cells.txt that build_kinker_input.py wrote.
# Usage: Rscript write_count_matrix.R INPUT_DIR
args<-commandArgs(trailingOnly=TRUE)
d<-args[1]
m<-Matrix::readMM(file.path(d,"counts.mtx"))
rownames(m)<-readLines(file.path(d,"genes.txt"))
colnames(m)<-readLines(file.path(d,"cells.txt"))
write.table(as.matrix(m),file=file.path(d,"count_matrix.txt"),sep="\t",quote=FALSE,row.names=TRUE)
cat("count_matrix.txt", nrow(m), "genes x", ncol(m), "cells\n")
