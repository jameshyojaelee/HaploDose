"""Compare tools on the Kinker single-line design from their evaluate.py line scores.

Each argument NAME=LINE_SCORES:METHOD takes rho_<METHOD> and sign_<METHOD> from one evaluate.py line_scores.tsv
(e.g. haplodose=.../line_scores.tsv:haplodose_log_ratio_gene). Lines present for any tool are listed
(unit = query line; specimen and line are the same in this design). The first argument is the reference:
paired Wilcoxon of its rho against each other tool, Holm over those pairs; effect = median paired difference
with a 95% percentile bootstrap CI over lines (10,000 resamples, seed 20261103). Lines missing for one tool are
dropped from that pair only (count reported as n_missing). Arms differ slightly between tools
(each run's used genes define the arm medians), so this compares each tool's profile against the same DepMap
truth on the arms it reports.
Usage: compare_tools.py DEST NAME=LINE_SCORES:METHOD [NAME=LINE_SCORES:METHOD ...]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


def main(dest, *specs):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    tables = {}
    for spec in specs:
        name, rest = spec.split('=', 1)
        path, method = rest.rsplit(':', 1)
        t = pd.read_csv(path, sep='\t').set_index('specimen')
        tables[name] = t[[f'rho_{method}', f'sign_{method}']].set_axis(['rho', 'sign'], axis=1)
    common = sorted(set.union(*(set(t.index) for t in tables.values())))
    rho = pd.DataFrame({n: t['rho'].reindex(common) for n, t in tables.items()})
    sign = pd.DataFrame({n: t['sign'].reindex(common) for n, t in tables.items()})
    rng = np.random.default_rng(20261103)
    rho.to_csv(dest / 'rho_by_line.tsv', sep='\t')
    ref = rho.columns[0]
    tests = []
    for other in rho.columns[1:]:
        d = (rho[ref] - rho[other]).dropna().to_numpy()
        boot = np.median(d[rng.integers(0, len(d), (10_000, len(d)))], axis=1)
        tests.append(dict(a=ref, b=other, n=len(d), n_missing=int(len(common) - len(d)),
                          median_a=float(rho[ref].median()), median_b=float(rho[other].median()),
                          median_diff=float(np.median(d)), ci95=[float(x) for x in np.quantile(boot, [0.025, 0.975])],
                          a_higher=int((d > 0).sum()), b_higher=int((d < 0).sum()), p=float(wilcoxon(d).pvalue)))
    if tests:
        holm = multipletests([t['p'] for t in tests], method='holm')[1]
        for t, p in zip(tests, holm):
            t['p_holm'] = float(p)
    summary = dict(n_lines=len(common), median_rho=rho.median().to_dict(),
                   mean_sign_informative=sign.mean().to_dict(), paired_tests=tests,
                   correction=f'Holm over the {len(tests)} pairs with {ref}; unit = cell line',
                   sources=list(specs))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
