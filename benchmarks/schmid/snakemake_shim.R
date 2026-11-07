# Minimal stand-in for the `snakemake` S4 object that Snakemake's script: directive injects, so the Schmid
# scripts (and the HaploDose-patched scorers) run with Rscript when snakemake is not installed.
# Usage, from their snakemake_pipeline/ directory:
#   source("path/to/snakemake_shim.R")
#   snakemake <- make_snakemake(input=list(...), output=list(...), params=list(...), log=list("x.log"))
#   source("scripts/evaluate_results_haplodose.R")
setClass("SnakemakeShim", representation(input="list", output="list", params="list", wildcards="list",
                                         threads="numeric", log="list", resources="list", config="list",
                                         rule="character", scriptdir="character", source="function"))

make_snakemake <- function(input=list(), output=list(), params=list(), wildcards=list(), threads=1,
                           log=list(), resources=list(), config=list(), rule="shim") {
  new("SnakemakeShim", input=input, output=output, params=params, wildcards=wildcards, threads=threads,
      log=log, resources=resources, config=config, rule=rule, scriptdir=getwd(), source=function(...) NULL)
}
