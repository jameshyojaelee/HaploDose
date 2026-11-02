#!/usr/bin/env bash
#SBATCH --job-name=build_cnv_field
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/build_cnv_field_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/build_cnv_field_%j.err
# Comparator env cnv_field (user-approved 2026-11-01): R 4.4 + CopyKAT and SCEVAN at pinned commits.
set -euo pipefail
unset PYTHONPATH R_LIBS R_LIBS_USER
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=/gpfs/commons/home/jameslee/micromamba/envs/cnv_field
test -d $ENV || $MM create -y -p $ENV -c conda-forge -c bioconda r-base=4.4 r-remotes r-matrix r-data.table \
  r-paralleldist r-dlm r-gplots r-rcolorbrewer r-mixtools r-cluster r-mcmcpack r-transport \
  r-pheatmap r-forcats r-dplyr r-ggplot2 r-rtsne r-ape r-tidytree r-ggrepel r-optparse \
  bioconductor-fgsea bioconductor-scran bioconductor-ggtree
$MM install -y -p $ENV -c conda-forge -c bioconda c-compiler cxx-compiler fortran-compiler make r-rcurl r-rcpphnsw \
  bioconductor-sparsematrixstats bioconductor-delayedmatrixstats r-doparallel r-foreach
export R_LIBS_USER=$ENV/lib/R/library PATH=$ENV/bin:$PATH
$ENV/bin/Rscript -e '
options(Ncpus = 8, repos = c(CRAN = "https://cloud.r-project.org"))
u <- function(repo, sha, ...) remotes::install_url(sprintf("https://github.com/%s/archive/%s.tar.gz", repo, sha), upgrade = "never", build = FALSE, ...)
u("miccec/yaGST", "56227df3ae183070c9d156af11c306ee799435e6", dependencies = FALSE)
u("navinlabcode/copykat", "01924c0999bb660064da8d3e1561615ac6efb76f")
u("AntonioDeFalco/SCEVAN", "5a49b88ac9445eeffcebb95404e3190992faac04")
library(copykat); library(SCEVAN)
writeLines(capture.output(sessionInfo()), "/gpfs/commons/home/jameslee/Haplodose/logs/cnv_field_sessionInfo.txt")
cat("OK\n")'
