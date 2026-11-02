"""Prespecified held-out analysis (docs/heldout_protocol.md), written and hashed before any held-out output existed.

E1 (primary): per-line Spearman rho vs DepMap arm truth from each tool's evaluate.py line_scores.tsv; Haplodose
column log_ratio_gene. Haplodose vs infercnvpy, R inferCNV, CopyKAT, SCEVAN and pseudobulk: paired two-sided
Wilcoxon over lines, Holm over these 5; median paired difference with a 95% percentile bootstrap CI (10,000 line
resamples, seed 20261103; ignores pool clustering, amendment 3a); lines without a profile from either tool are
dropped from that pair and counted per tool (amendment 3b). Superiority: Holm p < 0.05 and CI above 0.
Sensitivity (same tests): lines not in Kinker 2020. Descriptive: pool-rotation medians (6 reference sets) with a
Wilcoxon over them, Haplodose log_ratio (bins), sign agreement on informative arms, tumour clones per specimen,
coverage (query cells with a clone label / query cells sent).
E2 (secondary): implant design, ari_all per line (score_implant.py); Haplodose vs each of the 4 comparators,
paired Wilcoxon, Holm over 4, bootstrap CI; by floor stratum (implanted cells < 30 vs >= 30) and by setting,
descriptive.

Amendment 4 (2026-11-02, before any method run): failures follow the protocol instead of stopping the analysis.
A tool with no scored output keeps its place in each Holm family: E1 lines it did not return are NaN (dropped from
that pair, counted); a tool with no E2 output scores ari_all 0 on every line; a comparison with no testable pairs
enters Holm with p = 1 and is flagged. Tool failure logs (<tool>_failures.txt) are copied into the summary.
Usage: analyse.py ROOT INPUT_DIR
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

SEED = 20261103
COMPARATORS = ['infercnvpy', 'infercnv_r', 'copykat', 'scevan']
E1_FAMILY = ['infercnvpy', 'infercnv_r', 'copykat', 'scevan', 'pseudobulk']
SETTING = ['size', 'direction', 'fraction']


def wilcoxon_p(d):
    d = np.asarray(d, float)
    if len(d) == 0 or not np.any(d != 0):              # no pairs, or every difference 0: nothing to test
        return float('nan')
    return float(wilcoxon(d).pvalue)


def paired(a, b, rng):
    both = pd.concat([a, b], axis=1, keys=['a', 'b'])
    missing = dict(a=int(both.a.isna().sum()), b=int(both.b.isna().sum()))
    d = (both.a - both.b).dropna().to_numpy()
    if len(d):
        boot = np.median(d[rng.integers(0, len(d), (10_000, len(d)))], axis=1)
        ci, med = [float(x) for x in np.quantile(boot, [0.025, 0.975])], float(np.median(d))
    else:
        ci, med = [float('nan'), float('nan')], float('nan')
    return dict(n=len(d), missing_a=missing['a'], missing_b=missing['b'], median_diff=med, mean_diff=float(np.mean(d)) if len(d) else float('nan'),
                ci95=ci, a_higher=int((d > 0).sum()), b_higher=int((d < 0).sum()), p=wilcoxon_p(d))


def family(tests):
    """Holm over the full, fixed family; an untestable comparison enters with p = 1."""
    ps = [1.0 if np.isnan(t['p']) else t['p'] for t in tests]
    holm = multipletests(ps, method='holm')[1]
    for t, p in zip(tests, holm):
        t['testable'] = not np.isnan(t['p'])
        t['p_holm'] = float(p)
        t['superior'] = bool(t['testable'] and p < 0.05 and t['ci95'][0] > 0)
    return tests


def eval_dir(root, tool):
    return root / 'single' / f'eval_{tool}'


def line_scores(root, tool):
    f = eval_dir(root, tool) / 'line_scores.tsv'
    return pd.read_csv(f, sep='\t').set_index('specimen') if f.exists() else None


def e1(root, lines, rng):
    hd = line_scores(root, 'haplodose')
    if hd is None:
        raise SystemExit('primary endpoint not computable: no Haplodose line scores (see single/haplodose_failures.txt)')
    rho = {'haplodose': hd['rho_haplodose_log_ratio_gene'], 'pseudobulk': hd['rho_pseudobulk']}
    sign = {'haplodose': hd['sign_haplodose_log_ratio_gene'], 'pseudobulk': hd['sign_pseudobulk']}
    for t in COMPARATORS:
        s = line_scores(root, t)
        empty = pd.Series(dtype=float)
        rho[t] = s[f'rho_{t}_log_ratio'] if s is not None else empty
        sign[t] = s[f'sign_{t}_log_ratio'] if s is not None else empty
    rho, sign = pd.DataFrame(rho).reindex(lines), pd.DataFrame(sign).reindex(lines)
    tests = family([dict(a='haplodose', b=t, **paired(rho['haplodose'], rho[t], rng)) for t in E1_FAMILY])
    return rho, sign, hd, tests


def query_cells(inp):
    c = pd.concat([pd.read_csv(inp / f'cells_r{r}.tsv', sep='\t', usecols=['cell_id', 'design', 'specimen'], dtype=str)
                   for r in (0, 1)])
    return c[c.design.eq('single_line')].drop_duplicates('cell_id').set_index('cell_id').specimen


def coverage(root, tool, sent):
    """Per line: query cells with a clone label / query cells sent; a line with no output scores 0."""
    single = root / 'single'
    runs = ([single / f'haplodose_r{r}' for r in (0, 1)] if tool == 'haplodose'
            else sorted(single.glob(f'r*_pool*/{tool}_2*')))
    got = []
    for run in runs:
        if (run / 'qc.json').exists():
            c = pd.read_csv(run / 'cells.tsv.gz', sep='\t', usecols=['cell_id', 'design', 'clone'], dtype=str)
            c = c[c.design.eq('single_line') & c.cell_id.isin(sent.index)]
            got.append(c.set_index('cell_id').clone.notna() & c.set_index('cell_id').clone.ne('NA'))
    lab = pd.concat(got) if got else pd.Series(dtype=bool)
    lab = lab[~lab.index.duplicated()].reindex(sent.index).fillna(False).astype(bool)
    return lab.groupby(sent).mean()


def failures(root):
    out = {}
    for f in sorted(root.glob('*/*_failures.txt')):
        out[f'{f.parent.name}/{f.name}'] = f.read_text().splitlines()
    return out


def main(root, inp):
    root, inp = Path(root), Path(inp)
    out = root / 'analysis'
    out.mkdir(exist_ok=False)
    meta = pd.read_csv(inp / 'lines.tsv', sep='\t')
    chosen = meta[meta.chosen]
    all_lines = sorted(chosen.line)
    new_lines = sorted(chosen.loc[~chosen.in_kinker, 'line'])
    rng = np.random.default_rng(SEED)
    rho, sign, hd, tests = e1(root, all_lines, rng)
    _, _, _, tests_new = e1(root, new_lines, np.random.default_rng(SEED))
    rho.to_csv(out / 'e1_rho_by_line.tsv', sep='\t')
    rot = {}
    for r in (0, 1):                                   # each line is a query in exactly one rotation
        c = pd.read_csv(inp / f'cells_r{r}.tsv', sep='\t', usecols=['line', 'pool', 'design'], dtype=str)
        for line, p in c[c.design.eq('single_line')].drop_duplicates('line')[['line', 'pool']].itertuples(index=False):
            rot[line] = f'{p}_r{r}'
    pool_medians = rho.groupby(pd.Series(rot).reindex(all_lines).to_numpy()).median()
    pool_tests = [dict(a='haplodose', b=t, n_sets=int((pool_medians['haplodose'] - pool_medians[t]).notna().sum()),
                       p_nominal=wilcoxon_p((pool_medians['haplodose'] - pool_medians[t]).dropna()))
                  for t in E1_FAMILY]
    # E3, descriptive: tumour clones per single-line specimen (1 expected) and coverage.
    sent = query_cells(inp)
    cov = pd.DataFrame({t: coverage(root, t, sent) for t in ['haplodose'] + COMPARATORS}).reindex(all_lines)
    cov.to_csv(out / 'e3_coverage_by_line.tsv', sep='\t')
    clones = {}
    for t in ['haplodose'] + COMPARATORS:
        f = eval_dir(root, t) / 'specimens.tsv'
        if f.exists():
            s = pd.read_csv(f, sep='\t').set_index('specimen').reindex(all_lines)
            clones[t] = dict(n_with_output=int(s.tumour_clones.notna().sum()), median=float(s.tumour_clones.median()),
                             share_one_clone=float(s.tumour_clones.eq(1).sum() / len(all_lines)),
                             distribution=s.tumour_clones.value_counts().sort_index().to_dict())
        else:
            clones[t] = 'no output'
    summary = dict(
        e1=dict(n_lines=len(all_lines), median_rho=rho.median().to_dict(), mean_sign=sign.mean().to_dict(),
                haplodose_bins_median_rho=float(hd['rho_haplodose_log_ratio'].reindex(all_lines).median()),
                tests=tests, correction='Holm over 5 comparisons with Haplodose; unit = cell line; '
                'bootstrap CI over lines ignores pool clustering'),
        e1_sensitivity_not_in_kinker=dict(n_lines=len(new_lines), median_rho=rho.reindex(new_lines).median().to_dict(),
                                          tests=tests_new),
        e1_pool_medians=pool_medians.round(4).to_dict(orient='index'),
        e1_pool_tests=dict(tests=pool_tests, note='descriptive: Wilcoxon over the 6 pool-rotation reference sets, '
                           'nominal p, no correction; the smallest attainable two-sided p with 6 sets is 0.03125'),
        e3=dict(tumour_clones=clones, mean_coverage=cov.mean().to_dict(),
                lines_without_output={t: int(cov[t].eq(0).sum()) for t in cov}),
        failures=failures(root))
    imp_score = root / 'implant' / 'score'
    if (imp_score / 'per_specimen.tsv').exists():
        per = pd.read_csv(imp_score / 'per_specimen.tsv', sep='\t')
        ari = per.pivot_table(index='specimen', columns='tool', values='ari_all')
        absent = [t for t in ['haplodose'] + COMPARATORS if t not in ari]
        for t in absent:                                # no output at all: every specimen scores 0 (protocol)
            ari[t] = 0.0
        e2 = family([dict(a='haplodose', b=t, **paired(ari['haplodose'], ari[t], np.random.default_rng(SEED)))
                     for t in COMPARATORS])
        floor = per.groupby(['tool', 'below_floor']).ari_all.agg(['size', 'mean']).round(4).reset_index()
        setting = (per.groupby(['tool'] + SETTING).agg(n=('ari_all', 'size'), mean_ari_all=('ari_all', 'mean'),
                                                       region_found=('region_found', 'mean'))
                   .round(4).reset_index())
        setting.to_csv(out / 'e2_by_setting.tsv', sep='\t', index=False)
        summary['e2'] = dict(n_lines=int(ari.shape[0]), mean_ari_all=ari.mean().to_dict(), tools_without_output=absent,
                             tests=e2, by_floor=floor.to_dict(orient='records'),
                             by_setting=setting.to_dict(orient='records'),
                             correction='Holm over 4 comparisons with Haplodose; unit = cell line')
    else:
        summary['e2'] = f'not scored: {imp_score}/per_specimen.tsv missing'
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == '__main__':
    main(*sys.argv[1:])
