"""Write the HaploDose-patched copies of the Schmid scorers next to this file, plus a unified diff against the pinned
originals (external/benchmark_scrnaseq_cnv_callers, commit in external/README.md). Each edit is a literal
replacement that must match exactly once, so a changed upstream file fails loudly instead of patching silently.
Usage: python make_patched_scorers.py"""
import difflib
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / 'external/benchmark_scrnaseq_cnv_callers/snakemake_pipeline/scripts'

SOURCE_OLD = 'source("scripts/lib.R")\n'
SOURCE_NEW = '''source("scripts/lib.R")

#HaploDose reader (HaploDose repo: benchmarks/schmid/lib_haplodose.R)
haplodose_lib<-snakemake@params$haplodose_lib
if(is.null(haplodose_lib)){
  haplodose_lib<-"scripts/lib_haplodose.R"
}
source(haplodose_lib)
'''
INPUT_OLD = 'input_numbat_obj<-snakemake@input$numbat_obj\n'
INPUT_NEW = '''input_numbat_obj<-snakemake@input$numbat_obj

input_haplodose_genes<-snakemake@input$haplodose_genes
input_haplodose_cells<-snakemake@input$haplodose_cells
'''

COMBINED_BLOCK_OLD = '''#Save the combined result data frame for later processing
write.table(combined_range,file=output_merged_results,'''
COMBINED_BLOCK_NEW = '''#Check if HaploDose output files are provided
if(! is.null(input_haplodose_genes)){
  print("Load HaploDose results (clone profiles, pseudobulk over non-reference cells)")
  haplodose_results<-read_haplodose(input_haplodose_genes,input_haplodose_cells,
                                    input_annot,input_ref_groups)

  num_filtered_res<-rbind(num_filtered_res,
                          data.frame(method="HaploDose",
                                     annotated_genes=length(haplodose_results[[1]]),
                                     annotated_cells=haplodose_results[[2]]))

  combined_range<-combine_range_objects(combined_range,haplodose_results[[1]],
                                        method_colname="haplodose")
} else {
  print("Skipping HaploDose (no output file provided).")
}

#Save the combined result data frame for later processing
write.table(combined_range,file=output_merged_results,'''
NAMES_OLD = '''                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("wgs_mean","hb_cnv","honeybadger",
                         "infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_cnv","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))'''
NAMES_NEW = '''                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat", "HaploDose"),
                       c("wgs_mean","hb_cnv","honeybadger",
                         "infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_cnv","scevan",
                         "numbat","numbat_cnv", "CONICSmat", "haplodose"))'''


def permethod_block(src):
    """Insert the HaploDose block before the final write of res_per_method."""
    anchor = 'write.table(res_per_method,file=output_res'
    block = '''#Check if HaploDose output files are provided
if(! is.null(input_haplodose_genes)){
  print("Load HaploDose results")
  haplodose_results<-read_haplodose(input_haplodose_genes,input_haplodose_cells,
                                    input_annot,input_ref_groups)
  combined_range<-combine_range_objects(wgs_results,haplodose_results[[1]],
                                        method_colname="haplodose")
  res<-get_all_metrics(combined_range,"haplodose")
  res_per_method<-rbind(res_per_method,res)
} else {
  print("Skipping HaploDose (no output file provided).")
}

'''
    return [(anchor, block + anchor)]


EDITS = {
    'evaluate_results.R': [(SOURCE_OLD, SOURCE_NEW), (INPUT_OLD, INPUT_NEW), (COMBINED_BLOCK_OLD, COMBINED_BLOCK_NEW),
                           (NAMES_OLD, NAMES_NEW)],
    'evaluate_results_permethod.R': [(SOURCE_OLD, SOURCE_NEW), (INPUT_OLD, INPUT_NEW)] + permethod_block(None),
}


def main():
    diffs = []
    for name, edits in EDITS.items():
        old = (SRC / name).read_text()
        new = old
        for a, b in edits:
            if new.count(a) != 1:
                raise ValueError(f'{name}: anchor found {new.count(a)} times: {a[:60]!r}')
            new = new.replace(a, b)
        out = HERE / name.replace('.R', '_haplodose.R')
        out.write_text(new)
        diffs += difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                      f'a/snakemake_pipeline/scripts/{name}', f'b/snakemake_pipeline/scripts/{name}')
    (HERE / 'evaluate_results_haplodose.patch').write_text(''.join(diffs))


if __name__ == '__main__':
    main()
