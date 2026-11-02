"""Score bulk RNA arm-level copy number against U-CAN WGS FACETS, per patient and per arm.

Each method gives a table sample, arm, value. DNA per arm is log2(TCN / ploidy) (dna_truth.facets_arms).
Per tumour: Spearman rho across arms. Methods are compared per patient (the unit) with paired Wilcoxon
signed-rank tests, Holm-corrected over all method pairs. Per arm, across patients: Spearman rho with DNA and
ROC AUC for DNA gain and for DNA loss (FACETS TCN - ploidy >= 0.5 / <= -0.5).
Usage: score.py OUTDIR NAME=ARMS_TSV[:COLUMN] ...   (default column: value)
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ucan_spatial'))
from dna_truth import CEN, facets_arms  # noqa: E402


def load(spec):
    name, rest = spec.split('=', 1)
    path, col = (rest.rsplit(':', 1) + ['value'])[:2] if ':' in rest else (rest, 'value')
    t = pd.read_csv(path, sep='\t')
    return name, t.rename(columns={col: name})[['sample', 'arm', name]]


def main(outdir, *specs):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    methods = dict(load(s) for s in specs)
    rna = None
    for name, t in methods.items():
        rna = t if rna is None else rna.merge(t, on=['sample', 'arm'], how='outer')
    rna = rna[rna['sample'].str.endswith('.T')]
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    rows, meta, skipped = [], [], []
    for s in sorted(rna['sample'].unique()):
        patient = s.split('.')[2]
        try:
            dna, _, purity, ploidy, _ = facets_arms(patient, cen)
        except (FileNotFoundError, ValueError) as e:     # no VCF, or purity/ploidy NA in the header
            skipped.append(dict(sample=s, reason=str(e)[:80]))
            continue
        rows.append(dna.reset_index().assign(sample=s))
        meta.append(dict(sample=s, purity=purity, ploidy=ploidy))
    dna = pd.concat(rows)
    meta = pd.DataFrame(meta)
    pd.DataFrame(skipped, columns=['sample', 'reason']).to_csv(out / 'skipped.tsv', sep='\t', index=False)
    table = rna.merge(dna, on=['sample', 'arm'], how='inner')
    table.to_csv(out / 'arm_table.tsv.gz', sep='\t', index=False)
    names = list(methods)

    per = table.groupby('sample').apply(
        lambda g: pd.Series({m: spearmanr(g[m], g.rel_cn, nan_policy='omit')[0] for m in names}), include_groups=False)
    per = per.join(meta.set_index('sample'))
    per.to_csv(out / 'per_sample_rho.tsv', sep='\t')

    tests = []
    for a, b in itertools.combinations(names, 2):
        d = per[[a, b]].dropna()
        tests.append(dict(a=a, b=b, n=len(d), median_a=d[a].median(), median_b=d[b].median(),
                          a_higher=int((d[a] > d[b]).sum()), p=wilcoxon(d[a], d[b]).pvalue))
    tests = pd.DataFrame(tests)
    tests['p_holm'] = multipletests(tests.p, method='holm')[1]
    tests.to_csv(out / 'paired_tests.tsv', sep='\t', index=False)

    arm_rows = []
    for arm, g in table.groupby('arm'):
        r = dict(arm=arm, n=len(g), n_gain=int(g.dna_state.eq('gain').sum()), n_loss=int(g.dna_state.eq('loss').sum()))
        for m in names:
            x = g[[m, 'rel_cn', 'dna_state']].dropna()
            r[f'rho_{m}'] = spearmanr(x[m], x.rel_cn)[0]
            for st, sign in (('gain', 1), ('loss', -1)):
                y = x.dna_state.eq(st)
                r[f'auc_{st}_{m}'] = roc_auc_score(y, sign * x[m]) if 10 <= y.sum() <= len(y) - 10 else np.nan
        arm_rows.append(r)
    arms = pd.DataFrame(arm_rows)
    arms.to_csv(out / 'per_arm.tsv', sep='\t', index=False)

    summary = dict(n_tumours=int(len(per)), n_skipped=len(skipped), methods=names,
                   median_rho={m: float(per[m].median()) for m in names},
                   median_rho_purity_below_0_4={m: float(per.loc[per.purity < 0.4, m].median()) for m in names},
                   n_purity_below_0_4=int((per.purity < 0.4).sum()),
                   median_arm_auc_gain={m: float(arms[f'auc_gain_{m}'].median()) for m in names},
                   median_arm_auc_loss={m: float(arms[f'auc_loss_{m}'].median()) for m in names},
                   paired_tests=tests.to_dict(orient='records'),
                   correction='Holm over all method pairs; unit = patient (one tumour each)')
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps({k: summary[k] for k in summary if k != 'paired_tests'}, indent=1, default=float))
    print(tests.round(4).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
