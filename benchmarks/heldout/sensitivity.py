"""Prespecified sensitivity analyses of protocol amendment 5 (docs/heldout_protocol.md), written before any held-out
output existed. Run after analyse.py; the primary endpoints, tests and frozen scripts are unchanged. Truth, scoring,
the paired test and the Holm family are imported from the frozen scripts (evaluate.py, score_implant.py,
analyse.py), not copied.

(a) E1 on a common gene universe. Per rotation r and pool p: the genes used (`genes.tsv.gz`) by every tool that
    produced output there (Haplodose's rotation run, each comparator's r<r>_pool<p> run with qc.json). Every
    tool is rescored with evaluate.score_run on DepMap CN restricted to that universe, so truth (same formula),
    pseudobulk and tool arm profiles use the same genes and arms are eligible at >= 5 common genes. Same tests
    as E1 (analyse.paired: paired two-sided Wilcoxon, median paired difference with a 10,000-resample line
    bootstrap, seed 20261103) and the same Holm family of 5; a tool without output there is NaN for those lines.
    Decision rule of amendment 5: SOTA claim 1 additionally needs, against each of infercnvpy, R inferCNV,
    CopyKAT and SCEVAN, a mean paired difference > 0 with its bootstrap 95% CI above 0 (10,000 line resamples,
    seed 20261103, one fresh generator per comparison).
(b) E3 coverage restricted to query cells present in each tool's own matrix: a cell counts as covered only if
    it has a clone label and is in the tool's matrix (Haplodose and infercnvpy: every cell sent; R inferCNV:
    native/expr_cells.txt; CopyKAT, SCEVAN: native/cells.txt). Per line (cells sent as denominator; a line with
    no output scores 0) and per tool n_in_matrix / n_sent.
(c) E2 repeated with cells absent from a tool's matrix left unlabelled instead of clustered from zero
    profiles: their clone is set to NA and score_implant.score_run is applied unchanged (ari_all: unassigned
    cells one extra label; a specimen not returned scores 0; a tool with no output scores 0 on every line).
    Same tests as E2 (paired, Holm over 4).
(d) Lines share 6 reference sets (3 compositions x 2 rotations). `superior` is reported as conditional on
    independence of lines. Descriptive, for each Haplodose comparison of E1, (a) and E2: a cluster sign-flip
    permutation over reference sets of the per-set median paired difference (statistic: mean of the set
    medians; all 2^k patterns, two-sided; with 6 sets p >= 1/32), and a linear mixed model of the per-line
    difference with a random intercept per reference set (statsmodels MixedLM, REML).
Usage: sensitivity.py ROOT INPUT_DIR   (writes ROOT/analysis/sensitivity/)
"""
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'kinker'))
import analyse as A  # noqa: E402
import evaluate as EV  # noqa: E402
import score_implant as SI  # noqa: E402

SEED, N_BOOT = A.SEED, 10_000
PUBLISHED = A.COMPARATORS                               # the four published tools of SOTA claim 1
CONDITIONAL = 'conditional on independence of lines (lines share 6 reference sets; amendment 5d)'


# Runs and cell presence ------------------------------------------------------------------------------------
def _complete(dirs):
    return [d for d in sorted(dirs) if (d / 'qc.json').exists()]


def haplodose_runs(design_dir):
    return {r: design_dir / f'haplodose_r{r}' for r in (0, 1) if (design_dir / f'haplodose_r{r}' / 'qc.json').exists()}


def comparator_runs(design_dir, tool):
    """{(rotation, pool): latest completed run} from r<r>_pool<p>/<tool>_2*/ (the run_tool.sbatch layout)."""
    out = {}
    for d in _complete(design_dir.glob(f'r*_pool*/{tool}_2*')):
        r, p = d.parent.name.split('_pool')
        out[(int(r[1:]), p)] = d                       # sorted: a later run of the same pool wins
    return out


def present_cells(run, tool):
    """Cell ids in the tool's own matrix, or None if the run does not record them."""
    if tool in ('haplodose', 'infercnvpy'):            # every cell sent is in the matrix
        return set(pd.read_csv(run / 'cells.tsv.gz', sep='\t', usecols=['cell_id'], dtype=str).cell_id)
    f = run / 'native' / ('expr_cells.txt' if tool == 'infercnv_r' else 'cells.txt')
    return set(f.read_text().split('\n')[:-1]) if f.exists() else None


def rotation_sets(inp):
    """line -> reference set '<pool>_r<r>' and (r, pool) -> its query lines, from the input cell tables."""
    line_set, pool_lines = {}, {}
    for r in (0, 1):
        c = pd.read_csv(inp / f'cells_r{r}.tsv', sep='\t', usecols=['specimen', 'pool', 'design'], dtype=str)
        c = c[c.design.eq('single_line')].drop_duplicates('specimen')
        for spec, p in c[['specimen', 'pool']].itertuples(index=False):
            line_set[spec] = f'{p}_r{r}'
            pool_lines.setdefault((r, p), set()).add(spec)
    return line_set, pool_lines


# Statistics ------------------------------------------------------------------------------------------------
def boot_mean(a, b, seed=SEED):
    d = (pd.concat([a, b], axis=1, keys=['a', 'b']).pipe(lambda x: x.a - x.b)).dropna().to_numpy()
    if not len(d):
        return dict(mean_diff=float('nan'), mean_ci95=[float('nan'), float('nan')])
    rng = np.random.default_rng(seed)
    boot = d[rng.integers(0, len(d), (N_BOOT, len(d)))].mean(axis=1)
    return dict(mean_diff=float(d.mean()), mean_ci95=[float(x) for x in np.quantile(boot, [0.025, 0.975])])


def cluster_signflip(d, sets):
    """Two-sided sign-flip permutation over reference sets of the per-set median paired difference."""
    med = pd.Series(d).dropna().groupby(pd.Series(sets)).median().dropna()
    k = len(med)
    if k == 0:
        return dict(n_sets=0, statistic=float('nan'), p=float('nan'), min_attainable_p=float('nan'))
    m = med.to_numpy()
    t = float(m.mean())
    flips = np.array(list(itertools.product([1.0, -1.0], repeat=k)))
    perm = (flips * m[None]).mean(axis=1)
    p = float(np.mean(np.abs(perm) >= abs(t) - 1e-12))
    return dict(n_sets=k, set_medians={s: float(v) for s, v in med.items()}, statistic=t, p=p,
                min_attainable_p=float(2 / 2 ** k) if k else float('nan'))


def mixed_model(d, sets):
    """Random intercept per reference set, REML. With few sets the set variance is often at its 0 boundary, where
    statsmodels' default gradient optimizers fail or return NaN standard errors (seen on Kinker). Several
    optimizers are tried; among fits with a finite standard error the highest REML log-likelihood is kept,
    preferring converged ones. status is 'ok' only for a converged fit with a finite standard error."""
    df = pd.DataFrame({'d': pd.Series(d), 'set': pd.Series(sets)}).dropna()
    base = dict(n_lines=len(df), n_sets=int(df.set.nunique()))
    if len(df) < 3 or df.set.nunique() < 2:
        return dict(base, status='not fitted (too few lines or sets)')
    model = MixedLM(df.d.to_numpy(), np.ones((len(df), 1)), groups=df.set.to_numpy())
    fits = []
    for method in ('lbfgs', 'powell', 'nm', 'bfgs'):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            try:
                fit = model.fit(reml=True, method=method)
                se = float(fit.bse_fe[0])
                ci = np.asarray(fit.conf_int())[0]
                res = dict(method=method, converged=bool(fit.converged), llf=float(fit.llf),
                           intercept=float(fit.fe_params[0]), se=se, ci95=[float(ci[0]), float(ci[1])],
                           p_wald=float(fit.pvalues[0]), set_variance=float(np.asarray(fit.cov_re)[0, 0]),
                           residual_variance=float(fit.scale))
            except Exception as err:                # singular fits etc.: try the next optimizer
                res = dict(method=method, converged=False, llf=-np.inf, se=float('nan'), error=repr(err)[:200])
        res['warnings'] = sorted({str(x.message)[:200] for x in w})
        fits.append(res)
    finite = [f for f in fits if np.isfinite(f['se']) and np.isfinite(f['llf'])]
    if not finite:
        return dict(base, status='failed: no optimizer gave a finite standard error', attempts=fits)
    best = max(finite, key=lambda f: (f['converged'], f['llf']))
    status = 'ok' if best['converged'] else 'not converged (best finite fit reported)'
    return dict(base, status=status, **best, other_methods={f['method']: dict(converged=f['converged'], llf=f['llf'],
                                                                                 intercept=f.get('intercept'))
                                                            for f in fits if f is not best})


def dependence(table, line_set, family):
    """(d) for each Haplodose comparison in `family` (columns of `table`, index = line)."""
    sets = {l: line_set.get(l) for l in table.index}
    out = []
    for t in family:
        d = (table['haplodose'] - table[t]).to_dict()
        out.append(dict(a='haplodose', b=t, cluster_signflip=cluster_signflip(d, sets), mixed_model=mixed_model(d, sets)))
    return out


def conditional(tests):
    for t in tests:
        t['superior_conditional_on'] = CONDITIONAL
    return tests


# (a) E1 on a common gene universe ----------------------------------------------------------------------------
def e1_common(root, inp, lines, pool_lines):
    single = root / 'single'
    hd = haplodose_runs(single)
    comp = {t: comparator_runs(single, t) for t in A.COMPARATORS}
    cen = pd.read_csv(EV.CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    cn = pd.read_csv(inp / 'depmap_cn.tsv.gz', sep='\t', index_col=0)
    cn = cn[~cn.index.duplicated()]
    X, genes, cells = EV.load_counts(inp / 'counts.h5')
    used = lambda run: set(pd.read_csv(run / 'genes.tsv.gz', sep='\t', usecols=['gene']).gene)
    rho, universe, failed = {}, [], []
    for (r, p), specs in sorted(pool_lines.items()):
        if r not in hd:                                 # no Haplodose run: neither Haplodose nor pseudobulk
            universe.append(dict(rotation=r, pool=p, tools='', n_genes=0, note='no Haplodose output'))
            continue
        runs = {'haplodose': hd[r], **{t: comp[t][(r, p)] for t in A.COMPARATORS if (r, p) in comp[t]}}
        common = set.intersection(*(used(run) for run in runs.values()))
        universe.append(dict(rotation=r, pool=p, tools=','.join(runs), n_genes=len(common),
                             n_genes_with_cn=int(cn.index.isin(common).sum())))
        cn_u = cn[cn.index.isin(common)]
        for tool, run in runs.items():
            try:
                s = EV.score_run(inp, run, X, genes, cells, cn_u, cen, tool)[0]
            except Exception as err:                    # as a failed evaluate.py: the tool's lines stay NaN
                failed.append(dict(rotation=r, pool=p, tool=tool, run=str(run), error=repr(err)[:300]))
                continue
            if not len(s):
                continue
            s = s[s.specimen.isin(specs)].set_index('specimen')
            if tool == 'haplodose':
                cols = {'haplodose': 'rho_haplodose_log_ratio_gene', 'pseudobulk': 'rho_pseudobulk'}
                for k, c in cols.items():
                    rho.setdefault(k, {}).update(s[c].to_dict())
                rho.setdefault('n_arms', {}).update(s['n_arms'].to_dict())
            else:
                rho.setdefault(tool, {}).update(s[f'rho_{tool}_log_ratio'].to_dict())
    cols = ['haplodose'] + A.E1_FAMILY
    table = pd.DataFrame({c: pd.Series(rho.get(c, {}), dtype=float) for c in cols + ['n_arms']}).reindex(lines)
    if table['haplodose'].isna().all():
        raise SystemExit('(a) not computable: no Haplodose line scores')
    rng = np.random.default_rng(SEED)
    tests = A.family([dict(a='haplodose', b=t, **A.paired(table['haplodose'], table[t], rng)) for t in A.E1_FAMILY])
    for t in tests:
        t.update(boot_mean(table['haplodose'], table[t['b']]))
        t['sensitivity_rule_met'] = bool(t['mean_diff'] > 0 and t['mean_ci95'][0] > 0)
    return table, pd.DataFrame(universe), conditional(tests), failed


# (b) E3 coverage restricted to the tool's own matrix -------------------------------------------------------------
def coverage_in_matrix(root, inp):
    single = root / 'single'
    sent = A.query_cells(inp)
    by_line, counts, flags = {}, [], {}
    for tool in ['haplodose'] + A.COMPARATORS:
        runs = list(haplodose_runs(single).values()) if tool == 'haplodose' else list(comparator_runs(single, tool).values())
        lab, pres, unknown = [], set(), []
        for run in runs:
            c = pd.read_csv(run / 'cells.tsv.gz', sep='\t', usecols=['cell_id', 'design', 'clone'], dtype=str)
            c = c[c.design.eq('single_line') & c.cell_id.isin(sent.index)].set_index('cell_id')
            lab.append(c.clone.notna() & c.clone.ne('NA'))
            p = present_cells(run, tool)
            if p is None:
                unknown.append(str(run))
            else:
                pres |= p
        if unknown:                                     # presence not recorded: restricted coverage not computable
            flags[tool] = f'no matrix cell list in {unknown}'
            by_line[tool] = pd.Series(np.nan, index=sorted(set(sent)))
            counts.append(dict(tool=tool, n_sent=len(sent), n_in_matrix=np.nan, fraction_in_matrix=np.nan))
            continue
        lab = pd.concat(lab) if lab else pd.Series(dtype=bool)
        lab = lab[~lab.index.duplicated()].astype(bool).reindex(sent.index, fill_value=False)
        inm = pd.Series(sent.index.isin(list(pres)), index=sent.index)
        by_line[tool] = (lab & inm).groupby(sent).mean()
        counts.append(dict(tool=tool, n_sent=len(sent), n_in_matrix=int(inm.sum()),
                           fraction_in_matrix=float(inm.mean()), n_labelled=int(lab.sum()),
                           n_labelled_not_in_matrix=int((lab & ~inm).sum()), n_runs=len(runs)))
    return pd.DataFrame(by_line), pd.DataFrame(counts), flags


# (c) E2 with absent cells unlabelled ------------------------------------------------------------------------------
def e2_masked(root, work):
    imp_dir = root / 'implant' / 'input'
    if not (imp_dir / 'implants.tsv').exists():
        return None, f'not scored: {imp_dir}/implants.tsv missing', None
    imp = pd.read_csv(imp_dir / 'implants.tsv', sep='\t').set_index('specimen')
    implanted = pd.concat([pd.read_csv(imp_dir / f'cells_r{r}.tsv', sep='\t',
                                       usecols=['cell_id', 'specimen', 'design', 'implanted']) for r in (0, 1)])
    implanted = implanted[implanted.design.eq('single_line')].set_index('cell_id').implanted
    design = root / 'implant'
    rows, masked = [], []
    for tool in ['haplodose'] + A.COMPARATORS:
        runs = list(haplodose_runs(design).values()) if tool == 'haplodose' else list(comparator_runs(design, tool).values())
        for i, run in enumerate(runs):
            c = pd.read_csv(run / 'cells.tsv.gz', sep='\t', dtype=str)
            p = present_cells(run, tool)
            absent = np.zeros(len(c), bool) if p is None else ~c.cell_id.isin(p).to_numpy()
            q = c.design.eq('single_line').to_numpy() & c.cell_id.isin(implanted.index).to_numpy()
            masked.append(dict(tool=tool, run=str(run), query_cells=int(q.sum()), absent_query_cells=int((absent & q).sum()),
                               relabelled_to_NA=int((absent & q & c.clone.notna().to_numpy() & c.clone.ne('NA').to_numpy()).sum()),
                               matrix_list=p is not None))
            tmp = work / f'{tool}_{i}'
            tmp.mkdir(parents=True)
            c.loc[absent, 'clone'] = 'NA'
            c.to_csv(tmp / 'cells.tsv.gz', sep='\t', index=False, compression='gzip')
            (tmp / 'differences.tsv').symlink_to((run / 'differences.tsv').resolve())
            rows += SI.score_run(tool, tmp, imp, implanted)
    s = pd.DataFrame(rows)
    if len(s):                                           # as score_implant.main: a specimen not returned scores 0
        grid = pd.MultiIndex.from_product([s.tool.unique(), imp.index], names=['tool', 'specimen'])
        missing = grid.difference(pd.MultiIndex.from_frame(s[['tool', 'specimen']]))
        if len(missing):
            s = pd.concat([s, pd.DataFrame({'tool': missing.get_level_values(0), 'specimen': missing.get_level_values(1),
                                            'n_scored': 0, 'coverage': 0.0, 'ari_all': 0.0, 'clones': 0,
                                            'region_found': False, 'status': 'missing'})], ignore_index=True)
    ari = (s.pivot_table(index='specimen', columns='tool', values='ari_all') if len(s)
           else pd.DataFrame(index=imp.index)).reindex(imp.index)
    absent_tools = [t for t in ['haplodose'] + A.COMPARATORS if t not in ari]
    for t in absent_tools:                               # as analyse.py: no output at all scores 0 on every line
        ari[t] = 0.0
    tests = conditional(A.family([dict(a='haplodose', b=t, **A.paired(ari['haplodose'], ari[t], np.random.default_rng(SEED)))
                                  for t in A.COMPARATORS]))
    return ari, dict(n_lines=int(ari.shape[0]), mean_ari_all=ari.mean().to_dict(), tools_without_output=absent_tools,
                     masked=masked, tests=tests, correction='Holm over 4 comparisons with Haplodose; unit = cell line'), s


# Main -----------------------------------------------------------------------------------------------------------
def main(root, inp):
    root, inp = Path(root), Path(inp)
    if not (root / 'analysis').is_dir():
        raise SystemExit(f'{root}/analysis missing: run analyse.py first (this script must not create it)')
    out = root / 'analysis' / 'sensitivity'
    out.mkdir(exist_ok=False)
    meta = pd.read_csv(inp / 'lines.tsv', sep='\t')
    lines = sorted(meta[meta.chosen].line)
    line_set, pool_lines = rotation_sets(inp)
    summary = dict(protocol='docs/heldout_protocol.md amendment 5', seed=SEED, n_boot=N_BOOT, n_lines=len(lines))

    # Primary E1 as analyse.py computes it, labelled conditional (5d); used by the claim-1 rule below.
    rho, _, _, tests = A.e1(root, lines, np.random.default_rng(SEED))
    summary['e1_primary'] = dict(tests=conditional(tests), dependence=dependence(rho, line_set, A.E1_FAMILY))

    table, universe, tests_a, failed_a = e1_common(root, inp, lines, pool_lines)
    table.to_csv(out / 'a_e1_common_universe_rho_by_line.tsv', sep='\t')
    universe.to_csv(out / 'a_gene_universe.tsv', sep='\t', index=False)
    summary['a_e1_common_universe'] = dict(
        median_rho=table.drop(columns='n_arms').median().to_dict(), tests=tests_a,
        dependence=dependence(table, line_set, A.E1_FAMILY), universe=universe.to_dict(orient='records'),
        scoring_failures=failed_a,
        correction='Holm over 5 comparisons with Haplodose; unit = cell line')

    primary = {t['b']: t['superior'] for t in summary['e1_primary']['tests']}
    sens = {t['b']: t['sensitivity_rule_met'] for t in tests_a}
    summary['claim1_decision'] = dict(
        primary_superior=dict((t, primary[t]) for t in PUBLISHED),
        common_universe_mean_diff_ci_above_0=dict((t, sens[t]) for t in PUBLISHED),
        claim_holds=bool(all(primary[t] for t in PUBLISHED) and all(sens[t] for t in PUBLISHED)),
        rule='primary E1 superiority over all four published tools AND, on the common gene universe, mean paired '
             'difference > 0 with bootstrap 95% CI (of the mean) above 0 against each of them', note=CONDITIONAL)

    cov, counts, flags = coverage_in_matrix(root, inp)
    cov = cov.reindex(lines)
    cov.to_csv(out / 'b_coverage_in_matrix_by_line.tsv', sep='\t')
    counts.to_csv(out / 'b_cells_in_matrix.tsv', sep='\t', index=False)
    summary['b_coverage_in_matrix'] = dict(mean_coverage=cov.mean().to_dict(), cells=counts.to_dict(orient='records'),
                                           lines_without_output={t: int(cov[t].eq(0).sum()) for t in cov},
                                           not_computable=flags)

    work = out / 'c_e2_masked_runs'
    ari, e2, per = e2_masked(root, work)
    if ari is None:
        summary['c_e2_absent_unlabelled'] = e2
    else:
        per.to_csv(out / 'c_e2_absent_unlabelled_per_specimen.tsv', sep='\t', index=False)
        e2['dependence'] = dependence(ari, line_set, A.COMPARATORS)
        e2['masked_runs'] = str(work)                    # masked cells.tsv.gz copies, kept for audit
        summary['c_e2_absent_unlabelled'] = e2
    imp_score = root / 'implant' / 'score' / 'per_specimen.tsv'
    if imp_score.exists():                               # (d) on the primary E2 scores
        prim = pd.read_csv(imp_score, sep='\t').pivot_table(index='specimen', columns='tool', values='ari_all')
        for t in ['haplodose'] + A.COMPARATORS:
            if t not in prim:
                prim[t] = 0.0
        summary['e2_primary_dependence'] = dependence(prim, line_set, A.COMPARATORS)
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == '__main__':
    main(*sys.argv[1:])
