"""Compare tools on the Kinker single-line design from their evaluate.py line scores.

Version 2 reports paired-population summaries and mean-effect bootstrap inference,
resampling whole reference pools when supplied. The historical version-1 protocol
below remains documentation of archived results, not the current inference rule.

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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hardening.statistics import SCORE_VERSION, adjust_family, paired_mean_comparison


def main(dest, *specs):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    tables = {}
    reference_clusters = None
    for spec in specs:
        name, rest = spec.split('=', 1)
        path, method = rest.rsplit(':', 1)
        t = pd.read_csv(path, sep='\t').set_index('specimen')
        if t.index.has_duplicates:
            raise ValueError(f'duplicate specimen IDs in {path}')
        if name in tables:
            raise ValueError(f'duplicate tool name: {name}')
        if not tables and 'pool' in t:
            reference_clusters = t['pool']
        tables[name] = t[[f'rho_{method}', f'sign_{method}']].set_axis(['rho', 'sign'], axis=1)
    if not tables:
        raise ValueError('at least one tool table is required')
    common = sorted(set.union(*(set(t.index) for t in tables.values())))
    rho = pd.DataFrame({n: t['rho'].reindex(common) for n, t in tables.items()})
    sign = pd.DataFrame({n: t['sign'].reindex(common) for n, t in tables.items()})
    rho.to_csv(dest / 'rho_by_line.tsv', sep='\t')
    ref = rho.columns[0]
    tests = []
    for other in rho.columns[1:]:
        record = paired_mean_comparison(rho[ref], rho[other], seed=20261103,
                                        family_size=len(rho.columns) - 1,
                                        clusters=reference_clusters.reindex(common).to_numpy()
                                        if reference_clusters is not None else None)
        tests.append(dict(a=ref, b=other, **record))
    tests = adjust_family(tests)
    summary = dict(score_version=SCORE_VERSION, evidence_role='development or retrospective regression',
                   n_lines=len(common), median_rho=rho.median().to_dict(),
                   mean_sign_informative=sign.mean().to_dict(), paired_tests=tests,
                   correction=f'Holm over {len(tests)} mean-effect pairs; simultaneous basic-bootstrap intervals',
                   replication='whole reference pools when supplied; otherwise assumes independent lines',
                   uncertainty='conditional on supplied profiles; few reference pools limit inference',
                   sources=list(specs))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
