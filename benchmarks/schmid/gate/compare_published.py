"""Gate comparison (STATUS item 5, amendment d): our comparator rerun vs the published per-dataset values.

PUBLISHED is a TSV the team transcribes from the paper's Source Data (Schmid et al. 2025, Nat Commun 16:8777) with
columns dataset, method, metric, value; method keys as in their method_names (infercnv_cnv, infercnv_expr, casper,
copykat, scevan_cnv, scevan, numbat, numbat_cnv, CONICSmat); metric in pearson, auc_gains, auc_losses. No published
number is typed into this script. EVAL_DIR is results/output_<DATASET>/evaluation of the gate run.
Writes gate_comparison.tsv (ours, published, difference) to OUT. The pass rule is whatever STATUS fixes; this script
only reports the gaps. Usage: compare_published.py EVAL_DIR DATASET SUFFIX PUBLISHED OUT   (SUFFIX _wes or _cnvkit)
"""
import sys
from pathlib import Path

import pandas as pd


def main(eval_dir, dataset, suffix, published, out):
    ev = Path(eval_dir)
    corr = pd.read_csv(ev / f'evaluation_cnv_prediction_corr{suffix}.tsv', sep='\t')
    corr = corr[corr.method1.isin(['WGS', 'WES', 'scWGS', 'GATK'])]
    names = {'InferCNV (CNV)': 'infercnv_cnv', 'InferCNV (Expr)': 'infercnv_expr', 'CaSpER': 'casper', 'copyKat': 'copykat',
             'SCEVAN (CNV)': 'scevan_cnv', 'SCEVAN (Expr)': 'scevan', 'Numbat (Expr)': 'numbat', 'Numbat (CNV)': 'numbat_cnv',
             'CONICSmat': 'CONICSmat'}
    ours = [pd.DataFrame({'method': corr.method2.map(names), 'metric': 'pearson', 'ours': corr.pearson})]
    auc = pd.read_csv(ev / f'evaluation_cnv_prediction_auc{suffix}.txt', sep='\t')
    for m in ('auc_gains', 'auc_losses'):
        ours.append(pd.DataFrame({'method': auc.method, 'metric': m, 'ours': auc[m]}))
    ours = pd.concat(ours).dropna(subset=['method'])
    pub = pd.read_csv(published, sep='\t')
    pub = pub[pub.dataset.eq(dataset)].rename(columns={'value': 'published'})
    t = ours.merge(pub[['method', 'metric', 'published']], on=['method', 'metric'], how='outer')
    t['difference'] = t.ours - t.published
    t.insert(0, 'dataset', dataset)
    Path(out).mkdir(parents=True, exist_ok=True)
    t.to_csv(Path(out) / 'gate_comparison.tsv', sep='\t', index=False)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:6])
