#!/usr/bin/env python3
"""One table: 5 tools as columns (Haplodose, infercnvpy, prototype, R inferCNV, Numbat), dataset groups as
row blocks, plus time, peak RSS and CPUs per tool and block.

Haplodose roots (peer runs): U-CAN ucan_spatial_20261101_0911, Moorman moorman_20261101_0911,
cell lines celllines_20261101_0916, YCC ycc_20261101_0918/<S>/{haplodose,vs_numbat}.
Comparators: the latest collected_*.tsv of ROOT (collect.py) plus ycc vs_numbat summaries.
Labelling rules: prototype held-out false rate is NA (0 by construction); comparator subclonal split
recall comes from cn_state runs (row says so); inferCNV malignant and clone rows are "tool + our rule".
Haplodose CPUs come from its sbatch header (--cpus-per-task); Haplodose time/RSS from qc.json.
Usage: five_tool_table.py ROOT  -> ROOT/five_tool_table_<stamp>.tsv and .md
"""
import datetime as dt
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path('/gpfs/commons/home/jameslee/Haplodose')
UCAN, MOOR, CELL, YCC = (H / 'results' / p for p in
                         ['ucan_spatial_20261101_0911', 'moorman_20261101_0911', 'celllines_20261101_0916', 'ycc_20261101_0918'])
SBATCH = {'U-CAN': 'ucan_spatial/run_sample.sbatch', 'Moorman': 'moorman/run.sbatch',
          'cell lines': 'cell_lines/run.sbatch', 'YCC': 'ycc/run.sbatch'}
TOOLS = ['haplodose', 'infercnvpy', 'prototype', 'infercnv_r', 'numbat']
UCAN_S = ['U0007', 'U0046', 'U0178', 'U0179', 'U0209', 'U0245', 'U0247', 'U0285', 'U0354', 'U0383', 'U0396',
          'U0438', 'U0440', 'U0469', 'U0991', 'U2761', 'U2955', 'U2964', 'U2972', 'U2975']
YCC_S = ['YCC061', 'YCC097', 'YCC155', 'YCC162', 'YCC230']


def jload(p):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else {}


def cpus(block):
    m = re.search(r'--cpus-per-task=(\d+)', (H / 'benchmarks' / SBATCH[block]).read_text())
    return int(m.group(1))


def holdout_rate(run):
    c = pd.read_csv(Path(run) / 'cells.tsv.gz', sep='\t', usecols=lambda x: x in ('role', 'ref_split', 'malignant_call'))
    h = c.role.eq('reference') & c.ref_split.eq('holdout')
    return float(c.loc[h, 'malignant_call'].astype(bool).mean()) if h.any() else np.nan


def frac(num, den):
    return f'{int(np.nansum(num))}/{int(np.nansum(den))}'


def med(x):
    x = pd.Series(x, dtype=float).dropna()
    return round(float(x.median()), 3) if len(x) else 'NA'


def main(root):
    root = Path(root)
    col = pd.read_csv(sorted(root.glob('collected_*.tsv'))[-1], sep='\t')
    rows = []

    def add(block, metric, vals):
        rows.append(dict(block=block, metric=metric, **{t: vals.get(t, '') for t in TOOLS}))

    def resources(block, datasets):
        v_t, v_r, v_c = {}, {}, {}
        for t in TOOLS[1:]:
            g = col[col.tool.eq(t) & col.dataset.isin(datasets)]
            if not len(g):
                continue
            v_t[t] = f'{g.wall_seconds.sum() / 60:.1f} min total, max {g.wall_seconds.max() / 60:.1f}'
            if g.maxrss_gib_getrusage.notna().any():
                v_r[t] = f'{g.maxrss_gib_getrusage.max():.2f}'
            else:   # no getrusage in qc.json: sacct step MaxRSS (30 s sampling, so a lower bound)
                unit = {'K': 1, 'M': 2**10, 'G': 2**20}
                kib = [float(x[:-1]) * unit[x[-1]] if x[-1:] in unit else 0.0 for x in g.maxrss_sacct.astype(str)]
                v_r[t] = f'{max(kib) / 2**20:.2f} (sacct)'
            v_c[t] = int(g.cpus.max())
        return v_t, v_r, v_c

    # U-CAN
    hs = {s: jload(UCAN / s / 'dna/summary.json') for s in UCAN_S}
    hq = {s: jload(UCAN / s / 'haplodose/qc.json') for s in UCAN_S}
    u = col[col.dataset.isin(UCAN_S)]
    blk = 'U-CAN (20 sections vs WGS)'
    v = {'haplodose': len([s for s in hs if hs[s]])} | {t: int(u.tool.eq(t).sum()) for t in TOOLS[1:4]}
    add(blk, 'sections scored', v)
    v = {'haplodose': med([hs[s].get('rho_log_ratio') for s in hs])} | {t: med(u[u.tool.eq(t)].rho_arms) for t in TOOLS[1:4]}
    add(blk, 'median arm rho (log ratio)', v)
    add(blk, 'median arm rho (Haplodose log dose)', {'haplodose': med([hs[s].get('rho_log_dose') for s in hs])})
    v = {'haplodose': frac([hs[s].get('altered_arms_same_direction', 0) for s in hs], [hs[s].get('dna_altered_arms', 0) for s in hs])}
    v |= {t: frac(u[u.tool.eq(t)].altered_same_dir, u[u.tool.eq(t)].altered) for t in TOOLS[1:4]}
    add(blk, 'DNA-altered arms called, same direction', v)
    v = {'haplodose': frac([hs[s].get('neutral_arms_called_altered', 0) for s in hs], [hs[s].get('neutral_arms', 0) for s in hs])}
    v |= {t: frac(u[u.tool.eq(t)].neutral_called, u[u.tool.eq(t)].neutral) for t in TOOLS[1:4]}
    add(blk, 'DNA-neutral arms called altered', v)
    add(blk, 'subclonal segments split (Haplodose: statistical test)',
        {'haplodose': frac([hs[s].get('subclonal_segments_split', 0) for s in hs], [hs[s].get('subclonal_segments', 0) for s in hs])})
    add(blk, 'subclonal segments split (comparators: cn_state runs)',
        {t: frac(u[u.tool.eq(t)].subclonal_split_cn_state_runs, u[u.tool.eq(t)].subclonal) for t in TOOLS[1:4]})
    v = {'haplodose': med([hs[s].get('subclonal_split_fraction_expected_by_chance') for s in hs])}
    v |= {t: med(u[u.tool.eq(t)].split_expected_by_chance) for t in TOOLS[1:4]}
    add(blk, 'split fraction expected by chance (median)', v)
    v = {'haplodose': med([hs[s].get('tumour_clones') for s in hs])} | {t: med(u[u.tool.eq(t)].tumour_clones) for t in TOOLS[1:4]}
    add(blk, 'median tumour clones (inferCNVs: tool + our rule)', v)
    v = {'haplodose': med([holdout_rate(UCAN / s / 'haplodose') for s in UCAN_S])}
    v |= {t: med(u[u.tool.eq(t)].holdout_false_malignant) for t in ('infercnvpy', 'infercnv_r')} | {'prototype': 'NA'}
    add(blk, 'median held-out reference false-malignant rate (inferCNVs: tool + our rule)', v)
    rt, rr, rc = resources(blk, UCAN_S)
    rt['haplodose'] = f"{sum(q.get('wall_seconds', 0) for q in hq.values()) / 60:.1f} min total, max {max(q.get('wall_seconds', 0) for q in hq.values()) / 60:.1f}"
    rr['haplodose'] = f"{max(q.get('max_rss_gib', 0) for q in hq.values()):.2f}"
    rc['haplodose'] = cpus('U-CAN')
    add(blk, 'time', rt); add(blk, 'peak RSS GiB', rr); add(blk, 'CPUs', rc)

    # Cell lines
    blk = 'cell lines GSE149224 (HCT116 vs DepMap)'
    he = jload(CELL / 'eval/summary.json')
    cl = col[col.dataset.eq('celllines')].set_index('tool')
    ev = {t: jload(Path(cl.loc[t, 'run']) / 'eval/summary.json') for t in TOOLS[1:4] if t in cl.index and isinstance(cl.loc[t, 'run'], str)}
    ev['haplodose'] = he
    g = lambda k, f=lambda x: x: {t: f(e.get(k)) if e.get(k) is not None else 'NA' for t, e in ev.items()}
    add(blk, 'median arm rho HCT116', {t: round(e.get('rho_arms_median', {}).get('HCT116', np.nan), 3) for t, e in ev.items()})
    add(blk, 'mixture line ARI (median)', g('mixture_ari_line_median', lambda x: round(x, 3)))
    add(blk, 'tumour clones per mixture (tool + our rule for inferCNVs)', g('mixture_tumour_clones', lambda x: ','.join(map(str, x))))
    add(blk, 'tumour clones per single-line specimen', g('single_line_tumour_clones', lambda x: ','.join(map(str, x))))
    add(blk, 'single-line ARI clone vs dose', g('single_line_ari_dose', lambda x: ','.join(f'{y:.3f}' for y in x)))
    add(blk, 'single-line ARI clone vs batch', g('single_line_ari_batch', lambda x: ','.join(f'{y:.3f}' for y in x)))
    hc = {'haplodose': med([holdout_rate(CELL / 'haplodose')])} | {t: med([cl.loc[t, 'holdout_false_malignant']]) for t in ('infercnvpy', 'infercnv_r') if t in cl.index}
    add(blk, 'held-out reference false-malignant rate', hc | {'prototype': 'NA'})
    rt, rr, rc = resources(blk, ['celllines'])
    q = jload(CELL / 'haplodose/qc.json')
    rt['haplodose'] = f"{q.get('wall_seconds', 0) / 60:.1f} min"; rr['haplodose'] = f"{q.get('max_rss_gib', 0):.2f}"; rc['haplodose'] = cpus('cell lines')
    add(blk, 'time', rt); add(blk, 'peak RSS GiB', rr); add(blk, 'CPUs', rc)

    # Moorman
    blk = 'Moorman vs MSK-IMPACT'
    m = col[col.dataset.eq('Moorman')].set_index('tool')
    ms = {t: jload(Path(m.loc[t, 'run']) / 'dna/summary.json') for t in TOOLS[1:4] if t in m.index and isinstance(m.loc[t, 'run'], str)}
    ms['haplodose'] = jload(MOOR / 'dna/summary.json')
    add(blk, 'specimen-panel pairs (exact match / patients)', {t: f"{s.get('n_pairs')} ({s.get('n_exact')}/{s.get('n_exact_patients')})" for t, s in ms.items()})
    add(blk, 'median arm rho, exact pairs', {t: round(s.get('median_rho_exact', np.nan), 3) for t, s in ms.items()})
    add(blk, 'altered arms same direction', {t: s.get('altered_same_direction') for t, s in ms.items()})
    add(blk, 'neutral arms called altered', {t: s.get('neutral_called_altered') for t, s in ms.items()})
    rt, rr, rc = resources(blk, ['Moorman'])
    q = jload(MOOR / 'haplodose/qc.json')
    rt['haplodose'] = f"{q.get('wall_seconds', 0) / 60:.1f} min"; rr['haplodose'] = f"{q.get('max_rss_gib', 0):.2f}"; rc['haplodose'] = cpus('Moorman')
    add(blk, 'time', rt); add(blk, 'peak RSS GiB', rr); add(blk, 'CPUs', rc)

    # YCC: agreement with native Numbat (no DNA truth)
    blk = 'YCC (agreement with native Numbat, 5 samples)'
    ys = {'haplodose': {s: jload(YCC / s / 'vs_numbat/summary.json') for s in YCC_S}}
    for t in TOOLS[1:]:
        ys[t] = {}
        for s in YCC_S:
            runs = sorted(glob.glob(str(root / s / f'{t}_2*' / 'vs_numbat/summary.json')))
            ys[t][s] = jload(runs[-1]) if runs else {}
    per = lambda k: {t: ' '.join('NA' if d.get(k) is None else f'{d[k]:.2f}' for d in v.values()) for t, v in ys.items()}
    add(blk, 'per sample order', {t: ' '.join(YCC_S) for t in TOOLS})
    add(blk, 'malignant kappa (inferCNVs: tool + our rule)', per('malignant_kappa'))
    add(blk, 'tool malignant / Numbat tumour cells', {t: ' '.join('NA' if not d else f"{d.get('tool_malignant')}/{d.get('numbat_tumour')}" for d in v.values()) for t, v in ys.items()})
    add(blk, 'clone ARI on cells malignant in both', per('clone_ari'))
    add(blk, 'arm rho vs Numbat lnFC', per('rho_arms'))
    add(blk, 'arm state agreement', per('arm_state_agreement'))
    add(blk, 'note', {'numbat': 'comparator rerun of Numbat; YCC097 rerun: no CNV after entropy filter (NA)'})
    rt, rr, rc = resources(blk, YCC_S)
    hq = [jload(YCC / s / 'haplodose/qc.json') for s in YCC_S]
    rt['haplodose'] = f"{sum(q.get('wall_seconds', 0) for q in hq) / 60:.1f} min total, max {max(q.get('wall_seconds', 0) for q in hq) / 60:.1f}"
    rr['haplodose'] = f"{max(q.get('max_rss_gib', 0) for q in hq):.2f}"; rc['haplodose'] = cpus('YCC')
    add(blk, 'time', rt); add(blk, 'peak RSS GiB', rr); add(blk, 'CPUs', rc)

    t = pd.DataFrame(rows).rename(columns={'infercnv_r': 'R inferCNV', 'haplodose': 'Haplodose', 'numbat': 'Numbat'})
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    t.to_csv(root / f'five_tool_table_{stamp}.tsv', sep='\t', index=False)
    md = ['| ' + ' | '.join(t.columns) + ' |', '|' + '---|' * len(t.columns)]
    md += ['| ' + ' | '.join(str(x) for x in r) + ' |' for r in t.itertuples(index=False)]
    (root / f'five_tool_table_{stamp}.md').write_text('\n'.join(md) + '\n')
    print(root / f'five_tool_table_{stamp}.tsv')


if __name__ == '__main__':
    main(sys.argv[1])
