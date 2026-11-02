"""Head-to-head on Kinker: merge evaluate.py line_scores from several tools (same lines) and test per line.
Unit = cell line (one query per line); paired Wilcoxon on arm-level Spearman rho, Holm over tool pairs.
The pseudobulk column is taken from the first eval dir.
Usage: compare.py DEST EVAL_DIR...
"""
import itertools
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

KEYS = ['specimen', 'line', 'pool']


def main(dest, *evals):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    merged, tools = None, []
    for i, e in enumerate(evals):
        t = pd.read_csv(Path(e) / 'line_scores.tsv', sep='\t')
        cols = [c for c in t.columns if c.endswith('_log_ratio') and c.startswith(('rho_', 'sign_'))]
        if i == 0:
            cols += ['rho_pseudobulk', 'sign_pseudobulk', 'n_cells', 'informative_arms']
        tools += [c[4:] for c in cols if c.startswith('rho_')]
        merged = t[KEYS + cols] if merged is None else merged.merge(t[KEYS + cols], on=KEYS, validate='1:1')
    merged.to_csv(dest / 'line_scores_merged.tsv', sep='\t', index=False)
    rows = []
    for a, b in itertools.combinations(tools, 2):
        x = merged[[f'rho_{a}', f'rho_{b}']].dropna()
        rows.append(dict(a=a, b=b, n=len(x), median_a=x.iloc[:, 0].median(), median_b=x.iloc[:, 1].median(),
                         a_higher=int((x.iloc[:, 0] > x.iloc[:, 1]).sum()),
                         mean_sign_a=merged[f'sign_{a}'].mean(), mean_sign_b=merged[f'sign_{b}'].mean(),
                         p=wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue))
    tests = pd.DataFrame(rows)
    tests['p_holm'] = multipletests(tests.p, method='holm')[1]
    tests.to_csv(dest / 'paired_tests.tsv', sep='\t', index=False)
    print(f'lines {len(merged)}; Holm over {len(tests)} tool pairs; unit = cell line')
    print(tests.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:])
