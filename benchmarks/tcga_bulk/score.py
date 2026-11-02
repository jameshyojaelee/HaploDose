"""Score bulk RNA arm-level copy number against TCGA ABSOLUTE, per cancer type and pooled.

Truth per arm: length-weighted log2(max(Modal_Total_CN, 0.25) / ploidy) over the ABSOLUTE segments
(TCGA_mastercalls.abs_segtabs.fixed.txt; purity and ploidy from abs_tables_JSedit, call status 'called' only).
Segments are hg19; arms are cut at the hg19 p-arm acen band end (UCSC cytoBand) and matched to the RNA arms
(hg38 bins, hg38 centromeres) by arm name, so no liftover (as benchmarks/moorman/dna_truth.py). DNA state:
gain if mean TCN - ploidy >= 0.5, loss if <= -0.5 (benchmarks/ucan_spatial/dna_truth.py GAIN).
Metrics as benchmarks/ucan_bulk/score.py: per tumour Spearman rho across arms; methods compared per patient
(unit = patient, one primary tumour each) with paired Wilcoxon signed-rank tests, Holm-corrected over all
method pairs, separately within each cancer type and within the pooled set. Per arm across patients:
Spearman rho with DNA and ROC AUC for DNA gain and for DNA loss (>= 10 cases and >= 10 non-cases).
Usage: score.py OUTDIR TYPE=ROOT_TYPE_DIR ... --methods NAME=RELPATH[:COLUMN] ...
  each method file is ROOT_TYPE_DIR/RELPATH, a table sample, arm, value (default column value).
"""
import argparse
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

CN = Path('/gpfs/commons/home/jameslee/CRC/data/external/tcga_pancan_cn')
CYTO_HG19 = Path('/gpfs/commons/home/jameslee/CRC/data/external/bian_cnv_annotations/cytoBand.txt.gz')
GAIN = 0.5


def hg19_p_arm_end():
    with gzip.open(CYTO_HG19, 'rt') as fh:
        c = pd.read_csv(fh, sep='\t', header=None, names=['chromosome', 'start', 'end', 'band', 'stain'])
    return c[c.stain.eq('acen') & c.band.str.startswith('p')].set_index('chromosome').end


def absolute_arms(samples):
    tab = pd.read_csv(CN / 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t')
    tab = tab[tab['call status'].eq('called') & tab.array.isin(samples)].set_index('array')
    seg = pd.read_csv(CN / 'TCGA_mastercalls.abs_segtabs.fixed.txt', sep='\t',
                      usecols=['Sample', 'Chromosome', 'Start', 'End', 'Modal_Total_CN'])
    seg = seg[seg.Sample.isin(tab.index) & seg.Modal_Total_CN.notna()]
    seg['chromosome'] = 'chr' + seg.Chromosome.astype(float).astype(int).astype(str)
    pend = hg19_p_arm_end()
    seg = seg[seg.chromosome.isin(pend.index)]
    parts = []
    for side in ('p', 'q'):
        c = seg.chromosome.map(pend)
        s = seg.assign(start=np.maximum(seg.Start, c) if side == 'q' else seg.Start,
                       end=np.minimum(seg.End, c) if side == 'p' else seg.End)
        parts.append(s[s.end > s.start].assign(arm=lambda x: x.chromosome + side))
    s = pd.concat(parts)
    s['w'] = s.end - s.start
    s['ploidy'] = s.Sample.map(tab.ploidy)
    s['rel'] = np.log2(np.maximum(s.Modal_Total_CN, 0.25) / s.ploidy)
    s['tw'] = s.Modal_Total_CN * s.w
    s['rw'] = s.rel * s.w
    a = s.groupby(['Sample', 'arm'])[['w', 'tw', 'rw']].sum()
    dna = pd.DataFrame({'rel_cn': a.rw / a.w, 'tcn': a.tw / a.w}).reset_index().rename(columns={'Sample': 'sample'})
    dna['ploidy'] = dna['sample'].map(tab.ploidy)
    dna['dna_state'] = np.select([dna.tcn - dna.ploidy >= GAIN, dna.tcn - dna.ploidy <= -GAIN], ['gain', 'loss'], 'neutral')
    meta = tab[['purity', 'ploidy', 'Genome doublings']].rename(columns={'Genome doublings': 'genome_doublings'})
    return dna, meta.rename_axis('sample').reset_index()


def evaluate(table, per, names):
    tests = []
    for a, b in itertools.combinations(names, 2):
        d = per[[a, b]].dropna()
        p = wilcoxon(d[a], d[b]).pvalue if len(d) >= 5 and (d[a] != d[b]).any() else np.nan
        tests.append(dict(a=a, b=b, n=len(d), median_a=d[a].median(), median_b=d[b].median(),
                          a_higher=int((d[a] > d[b]).sum()), p=p))
    tests = pd.DataFrame(tests)
    ok = tests.p.notna()
    tests['p_holm'] = np.nan
    if ok.any():
        tests.loc[ok, 'p_holm'] = multipletests(tests.p[ok], method='holm')[1]
    arm_rows = []
    for arm, g in table.groupby('arm'):
        r = dict(arm=arm, n=len(g), n_gain=int(g.dna_state.eq('gain').sum()), n_loss=int(g.dna_state.eq('loss').sum()))
        for m in names:
            x = g[[m, 'rel_cn', 'dna_state']].dropna()
            r[f'rho_{m}'] = spearmanr(x[m], x.rel_cn)[0] if len(x) > 3 else np.nan
            for st, sign in (('gain', 1), ('loss', -1)):
                y = x.dna_state.eq(st)
                r[f'auc_{st}_{m}'] = roc_auc_score(y, sign * x[m]) if 10 <= y.sum() <= len(y) - 10 else np.nan
        arm_rows.append(r)
    return tests, pd.DataFrame(arm_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('types', nargs='+', help='TYPE=DIR')
    ap.add_argument('--methods', nargs='+', required=True, help='NAME=RELPATH[:COLUMN]')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    methods = {}
    for spec in args.methods:
        name, rest = spec.split('=', 1)
        path, col = rest.rsplit(':', 1) if ':' in rest else (rest, 'value')
        methods[name] = (path, col)
    names = list(methods)

    rna = []
    for spec in args.types:
        t, d = spec.split('=', 1)
        s = pd.read_csv(Path(d) / 'samples.tsv', sep='\t')
        tum = set(s.loc[s.role.eq('tumour'), 'sample'])
        m = None
        for name, (path, col) in methods.items():
            x = pd.read_csv(Path(d) / path, sep='\t').rename(columns={col: name})[['sample', 'arm', name]]
            m = x if m is None else m.merge(x, on=['sample', 'arm'], how='outer')
        rna.append(m[m['sample'].isin(tum)].assign(type=t))
    rna = pd.concat(rna, ignore_index=True)
    dna, meta = absolute_arms(set(rna['sample']))
    table = rna.merge(dna, on=['sample', 'arm'], how='inner')
    table.to_csv(out / 'arm_table.tsv.gz', sep='\t', index=False)
    no_truth = sorted(set(rna['sample']) - set(table['sample']))
    pd.Series(no_truth, name='sample').to_csv(out / 'tumours_without_absolute.tsv', sep='\t', index=False)

    per = table.groupby(['type', 'sample']).apply(
        lambda g: pd.Series({m: spearmanr(g[m], g.rel_cn, nan_policy='omit')[0] for m in names}), include_groups=False)
    per = per.reset_index().merge(meta, on='sample', how='left')
    per.to_csv(out / 'per_sample_rho.tsv', sep='\t', index=False)

    all_tests, all_arms, summary = [], [], []
    groups = [(t, g) for t, g in table.groupby('type')] + [('pooled', table)]
    for t, g in groups:
        p = per if t == 'pooled' else per[per.type.eq(t)]
        tests, arms = evaluate(g, p.set_index('sample'), names)
        all_tests.append(tests.assign(type=t))
        all_arms.append(arms.assign(type=t))
        row = dict(type=t, n_tumours=len(p), n_purity_below_0_4=int((p.purity < 0.4).sum()))
        for m in names:
            row[f'median_rho_{m}'] = p[m].median()
            row[f'median_rho_purity_below_0_4_{m}'] = p.loc[p.purity < 0.4, m].median()
            row[f'median_arm_auc_gain_{m}'] = arms[f'auc_gain_{m}'].median()
            row[f'median_arm_auc_loss_{m}'] = arms[f'auc_loss_{m}'].median()
            row[f'arms_auc_gain_{m}'] = int(arms[f'auc_gain_{m}'].notna().sum())
        summary.append(row)
    pd.concat(all_tests).to_csv(out / 'paired_tests.tsv', sep='\t', index=False)
    pd.concat(all_arms).to_csv(out / 'per_arm.tsv', sep='\t', index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(out / 'summary.tsv', sep='\t', index=False)
    (out / 'summary.json').write_text(json.dumps(dict(
        methods=names, n_tumours_scored=int(len(per)), n_tumours_without_absolute=len(no_truth),
        correction='Holm over all method pairs, within each cancer type and within pooled; unit = patient',
        truth='ABSOLUTE log2(TCN/ploidy) per arm, hg19 arms matched by name'), indent=1))
    pd.set_option('display.width', 250)
    print(summary[['type', 'n_tumours'] + [f'median_rho_{m}' for m in names]].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
