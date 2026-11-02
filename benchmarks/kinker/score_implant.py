"""Score the Kinker subclone design (implant.py). Unit = query line (one implanted specimen per line).
Per tool and specimen, over the specimen's query cells:
- ARI between the tool's clone labels and the implanted / not implanted split;
- clones = number of distinct labels (2 expected);
- region found = some clone difference in differences.tsv overlaps the implanted region (same chromosome).
Paired Wilcoxon on ARI between tools, Holm over tool pairs; per-setting means are descriptive.
Usage: score_implant.py IMPLANT_DIR DEST TOOL=RUN_DIR [TOOL=RUN_DIR ...]  (a tool may repeat, one per run)
"""
import itertools
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

SETTING = ['size', 'direction', 'fraction']


def score_run(tool, run, imp, implanted):
    run = Path(run)
    cells = pd.read_csv(run / 'cells.tsv.gz', sep='\t', usecols=['cell_id', 'specimen', 'design', 'clone'], dtype=str)
    cells = cells[cells.design.eq('single_line') & cells.specimen.isin(imp.index)]
    diff = pd.read_csv(run / 'differences.tsv', sep='\t')
    rows = []
    for spec, c in cells.groupby('specimen'):
        c = c[c.clone.notna() & c.clone.ne('NA')]
        y = implanted.loc[c.cell_id].to_numpy()
        e = imp.loc[spec]
        d = diff[diff.specimen.eq(spec) & diff.chromosome.eq(e.chromosome)]
        found = bool(((d.start <= e.end) & (d.end >= e.start)).any())
        rows.append(dict(tool=tool, specimen=spec, n_scored=len(c), ari=adjusted_rand_score(y, c.clone),
                         clones=c.clone.nunique(), region_found=found))
    return rows


def main(imp_dir, dest, *runs):
    imp_dir, dest = Path(imp_dir), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    imp = pd.read_csv(imp_dir / 'implants.tsv', sep='\t').set_index('specimen')
    implanted = pd.concat([pd.read_csv(imp_dir / f'cells_r{r}.tsv', sep='\t', usecols=['cell_id', 'specimen', 'design', 'implanted'])
                           for r in (0, 1)])
    implanted = implanted[implanted.design.eq('single_line')].set_index('cell_id').implanted
    rows = []
    for arg in runs:
        tool, run = arg.split('=', 1)
        rows += score_run(tool, run, imp, implanted)
    s = pd.DataFrame(rows).merge(imp[SETTING + ['n_cells', 'n_implanted']], left_on='specimen', right_index=True)
    s.to_csv(dest / 'per_specimen.tsv', sep='\t', index=False)
    by = s.groupby(['tool'] + SETTING).agg(n=('ari', 'size'), mean_ari=('ari', 'mean'), median_clones=('clones', 'median'),
                                            region_found=('region_found', 'mean')).reset_index()
    by.to_csv(dest / 'by_setting.tsv', sep='\t', index=False)
    wide = s.pivot(index='specimen', columns='tool', values='ari')
    tests = []
    for a, b in itertools.combinations(wide.columns, 2):
        x = wide[[a, b]].dropna()
        tests.append(dict(a=a, b=b, n=len(x), mean_a=x[a].mean(), mean_b=x[b].mean(), a_higher=int((x[a] > x[b]).sum()),
                          b_higher=int((x[b] > x[a]).sum()), p=wilcoxon(x[a], x[b]).pvalue))
    tests = pd.DataFrame(tests)
    if len(tests):
        tests['p_holm'] = multipletests(tests.p, method='holm')[1]
    tests.to_csv(dest / 'paired_tests.tsv', sep='\t', index=False)
    summary = dict(n_specimens=int(wide.shape[0]),
                   overall=s.groupby('tool').agg(mean_ari=('ari', 'mean'), ari_ge_0_5=('ari', lambda a: float((a >= 0.5).mean())),
                                                 median_clones=('clones', 'median'),
                                                 region_found=('region_found', 'mean')).to_dict(orient='index'),
                   paired_tests=tests.to_dict(orient='records'),
                   correction='Holm over tool pairs; unit = cell line (one implanted specimen each)')
    (dest / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
    print(by.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:])
