"""Score the Kinker subclone design (implant.py). Unit = query line (one implanted specimen per line).
Per tool and specimen, over the specimen's query cells:
- ari: ARI between the tool's clone labels and the implanted / not implanted split, assigned cells only (NaN
  when fewer than 2 cells are assigned; sklearn gives 1.0 on empty input);
- ari_all (primary, failure-aware): the same over all query cells, unassigned cells as one extra label; a
  specimen missing from a tool's output scores 0. coverage = assigned / query cells;
- clones = number of distinct labels (2 expected);
- region found = some clone difference in differences.tsv overlaps the implanted region (same chromosome).
Strata (descriptive): setting, and below / at or above the proposal size floor (n_implanted < MIN_CELLS).
Paired Wilcoxon on ari_all between tools, Holm over tool pairs; per-setting means are descriptive.
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
# Proposal size floor: CloneParams.min_cells default. A split proposal needs at least this many cells per side, so
# an implanted minority below it cannot be proposed as is (review 2026-11-02 #5). Reported as a stratum; the
# primary endpoint keeps every specimen.
MIN_CELLS = 30


def score_run(tool, run, imp, implanted):
    run = Path(run)
    cells = pd.read_csv(run / 'cells.tsv.gz', sep='\t', usecols=['cell_id', 'specimen', 'design', 'clone'], dtype=str)
    cells = cells[cells.design.eq('single_line') & cells.specimen.isin(imp.index)]
    diff = pd.read_csv(run / 'differences.tsv', sep='\t')
    rows = []
    for spec, c_all in cells.groupby('specimen'):
        ok = c_all.clone.notna() & c_all.clone.ne('NA')
        c = c_all[ok]
        y = implanted.loc[c.cell_id].to_numpy()
        y_all = implanted.loc[c_all.cell_id].to_numpy()
        e = imp.loc[spec]
        d = diff[diff.specimen.eq(spec) & diff.chromosome.eq(e.chromosome)]
        found = bool(((d.start <= e.end) & (d.end >= e.start)).any())
        rows.append(dict(tool=tool, specimen=spec, n_query=len(c_all), n_scored=len(c), coverage=len(c) / len(c_all),
                         ari=adjusted_rand_score(y, c.clone) if len(c) >= 2 else float('nan'),
                         ari_all=adjusted_rand_score(y_all, c_all.clone.where(ok, '__unassigned__')),
                         clones=c.clone.nunique(), region_found=found, status='ok' if len(c) >= 2 else 'no_output'))
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
    s = pd.DataFrame(rows)
    # A specimen a tool did not return at all is a failure, not a drop-out: ari_all 0, coverage 0.
    grid = pd.MultiIndex.from_product([s.tool.unique(), imp.index], names=['tool', 'specimen'])
    missing = grid.difference(pd.MultiIndex.from_frame(s[['tool', 'specimen']]))
    if len(missing):
        s = pd.concat([s, pd.DataFrame({'tool': missing.get_level_values(0), 'specimen': missing.get_level_values(1),
                                        'n_scored': 0, 'coverage': 0.0, 'ari_all': 0.0, 'clones': 0,
                                        'region_found': False, 'status': 'missing'})], ignore_index=True)
    s = s.merge(imp[SETTING + ['n_cells', 'n_implanted']], left_on='specimen', right_index=True)
    s['below_floor'] = s.n_implanted < MIN_CELLS
    s.to_csv(dest / 'per_specimen.tsv', sep='\t', index=False)
    by = s.groupby(['tool'] + SETTING + ['below_floor']).agg(n=('ari_all', 'size'), mean_ari_all=('ari_all', 'mean'), mean_ari=('ari', 'mean'),
                                            median_clones=('clones', 'median'),
                                            region_found=('region_found', 'mean')).reset_index()
    by.to_csv(dest / 'by_setting.tsv', sep='\t', index=False)
    wide = s.pivot(index='specimen', columns='tool', values='ari_all')
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
                   overall=s.groupby('tool').agg(mean_ari_all=('ari_all', 'mean'), mean_ari=('ari', 'mean'),
                                                 ari_ge_0_5=('ari_all', lambda a: float((a >= 0.5).mean())),
                                                 mean_coverage=('coverage', 'mean'),
                                                 n_no_output=('status', lambda x: int(x.ne('ok').sum())),
                                                 median_clones=('clones', 'median'),
                                                 region_found=('region_found', 'mean')).to_dict(orient='index'),
                   by_floor={f'{t}|{"below" if b else "at_or_above"}_{MIN_CELLS}': dict(
                       n=int(len(g)), mean_ari_all=float(g.ari_all.mean()), region_found=float(g.region_found.mean()))
                       for (t, b), g in s.groupby(['tool', 'below_floor'])},
                   n_below_floor=int(s.drop_duplicates('specimen').below_floor.sum()),
                   paired_tests=tests.to_dict(orient='records'),
                   endpoint='ari_all (failure-aware: unassigned cells one extra label, missing specimen 0)',
                   correction='Holm over tool pairs; unit = cell line (one implanted specimen each)')
    (dest / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
    print(by.to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:])
