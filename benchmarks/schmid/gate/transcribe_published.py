"""Transcribe the published per-output gate values (claim 5 amendment (g), docs/STATUS.md) from the Schmid et al.
2025 Source Data (Nat Commun 16:8777, 41467_2025_62359_MOESM11_ESM.zip, md5 ef0ef583893661e34c3bcdf3eac3de21).

Source file: figure_source_data/suppfigure_4_metrics_extended.tsv (the zip holds TSVs, no sheets). Its
'AUC (gain)'/'AUC (loss)' are the ROCR full AUCs our scorer writes as auc_gains/auc_losses (lib.R:1309, 1325);
its norm_value is (AUC - 0.5) * 2 floored at 0, as in the paper's summary script. Figure 2A shows the pROC partial
AUC instead, so it is not used for the gate. 'value' (not 'norm_value') is transcribed. Correlation is rounded to
3 decimals in the source.
Method names -> keys as in evaluate_results.R:279-286.
Usage: python -I transcribe_published.py SOURCE_TSV OUT_TSV
"""
import csv
import sys

DATASETS = {'COLO320': 'primary', 'BCC06': 'primary', 'BCC06post': 'primary', 'MM': 'primary',
            'A375': 'line_overlap', 'MCF7': 'line_overlap', 'HCT116': 'line_overlap'}
METRICS = {'Correlation': 'pearson', 'AUC (gain)': 'auc_gains', 'AUC (loss)': 'auc_losses'}
METHODS = {'InferCNV (CNV)': 'infercnv_cnv', 'InferCNV (Expr)': 'infercnv_expr', 'CaSpER': 'casper',
           'copyKat': 'copykat', 'SCEVAN (CNV)': 'scevan_cnv', 'SCEVAN (Expr)': 'scevan',
           'Numbat (Expr)': 'numbat', 'Numbat (CNV)': 'numbat_cnv', 'CONICSmat': 'CONICSmat'}
SOURCE_FILE = '41467_2025_62359_MOESM11_ESM.zip:figure_source_data/suppfigure_4_metrics_extended.tsv'


def main(src, out):
    rows, seen = [], set()
    with open(src, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['dataset'] not in DATASETS or r['metric'] not in METRICS:
                continue
            key = (r['dataset'], METHODS[r['method']], METRICS[r['metric']])
            if key in seen:
                raise ValueError(f'duplicate row {key}')
            seen.add(key)
            rows.append(dict(dataset=r['dataset'], set=DATASETS[r['dataset']], method=key[1], metric=key[2],
                             value=r['value'], published_method=r['method'], published_metric=r['metric'],
                             published_norm_value=r['norm_value'], source_file=SOURCE_FILE, sheet='NA (tsv)'))
    missing = [(d, m) for d in DATASETS for m in METHODS.values()
               if not any(x['dataset'] == d and x['method'] == m for x in rows)]
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} rows; outputs absent from the source: {missing}')


if __name__ == '__main__':
    main(*sys.argv[1:3])
