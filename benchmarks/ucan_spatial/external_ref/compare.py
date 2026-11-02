"""Baseline (immune + stromal reference) vs external normal-epithelium reference, U-CAN Visium HD vs FACETS.

Usage: compare.py BASELINE_ROOT EXTREF_ROOT
Patient is the unit (20, one section each). One test only: paired Wilcoxon signed-rank on per-patient
Spearman rho (RNA arm log ratio vs DNA log2 TCN/ploidy), on patients with rho defined in both runs.
Writes EXTREF_ROOT/comparison_per_patient.tsv and comparison_totals.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

base, ext = Path(sys.argv[1]), Path(sys.argv[2])


def load(root):
    t = pd.read_csv(root / 'cohort_summary_dna.tsv', sep='\t').set_index('sample')
    spec, rc, arms = [], [], []
    for s in t.index:
        sp = pd.read_csv(root / s / 'haplodose/specimens.tsv', sep='\t')
        spec.append(dict(sample=s, tau=float(sp.tau.iloc[0]), mixture=str(sp.mixture.iloc[0])))
        r = pd.read_csv(root / s / 'haplodose/reference_check.tsv', sep='\t')
        rc.append(dict(sample=s, ref_tested=int(r.n_tested.sum()), ref_called=int(r.n_called_aneuploid.sum()),
                       ref_expected=float((r.expected_fraction * r.n_tested).sum())))
        arms.append(pd.read_csv(root / s / 'dna/arms.tsv', sep='\t'))
    t = t.join(pd.DataFrame(spec).set_index('sample')).join(pd.DataFrame(rc).set_index('sample'))
    t['ref_false_rate'] = t.ref_called / t.ref_tested
    t['subclonal_chance'] = t.subclonal_split_fraction_expected_by_chance.fillna(0) * t.subclonal_segments
    return t, pd.concat(arms, ignore_index=True)


b, ab = load(base)
e, ae = load(ext)
cols = ['tumour_clones', 'rho_log_ratio', 'altered_arms_same_direction', 'dna_altered_arms',
        'neutral_arms_called_altered', 'neutral_arms', 'subclonal_segments_split', 'subclonal_chance',
        'subclonal_segments', 'ref_false_rate', 'tau', 'mixture']
per = b[cols].add_suffix('_base').join(e[cols].add_suffix('_ext'))
per['rho_diff'] = per.rho_log_ratio_ext - per.rho_log_ratio_base
per.to_csv(ext / 'comparison_per_patient.tsv', sep='\t')

both = per.dropna(subset=['rho_log_ratio_base', 'rho_log_ratio_ext'])
w = wilcoxon(both.rho_log_ratio_ext, both.rho_log_ratio_base)
arms = ab[['sample', 'arm', 'rel_cn', 'rna_log_ratio']].merge(
    ae[['sample', 'arm', 'rna_log_ratio']], on=['sample', 'arm'], suffixes=('_base', '_ext'))
common = arms.dropna(subset=['rel_cn', 'rna_log_ratio_base', 'rna_log_ratio_ext'])


def pooled(a, col):
    a = a.dropna(subset=['rel_cn', col])
    return spearmanr(a[col], a.rel_cn).statistic, len(a)


rows = []
for name, t, a, col in [('baseline', b, ab, 'rna_log_ratio'), ('external_ref', e, ae, 'rna_log_ratio')]:
    r_all, n_all = pooled(a, col)
    r_com = spearmanr(common[f'{col}_{"base" if name == "baseline" else "ext"}'], common.rel_cn).statistic
    rows.append(dict(run=name, patients_with_tumour_clone=int((t.tumour_clones > 0).sum()),
                     median_rho=t.rho_log_ratio.median(), median_rho_18_common=t.loc[both.index].rho_log_ratio.median(),
                     pooled_rho_all_arms=r_all, pooled_n_arms=n_all, pooled_rho_common_arms=r_com,
                     n_common_arms=len(common),
                     altered_right=f'{int(t.altered_arms_same_direction.sum())}/{int(t.dna_altered_arms.sum())}',
                     neutral_called=f'{int(t.neutral_arms_called_altered.sum())}/{int(t.neutral_arms.sum())}',
                     tumour_clones_total=int(t.tumour_clones.sum()), tumour_clones_median=t.tumour_clones.median(),
                     subclonal_split=f'{int(t.subclonal_segments_split.sum())}/{int(t.subclonal_segments.sum())}',
                     subclonal_chance=round(t.subclonal_chance.sum(), 1),
                     ref_holdout_false=f'{int(t.ref_called.sum())}/{int(t.ref_tested.sum())}',
                     ref_holdout_expected=round(t.ref_expected.sum(), 1), median_tau=t.tau.median()))
tot = pd.DataFrame(rows)
tot.to_csv(ext / 'comparison_totals.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(per.drop(columns=['mixture_base', 'mixture_ext']).round(3).to_string())
print(tot.T.to_string())
print(f'\nPaired Wilcoxon signed-rank on per-patient rho (ext vs base), n = {len(both)} patients with rho in both '
      f'runs: statistic {w.statistic:.1f}, two-sided p = {w.pvalue:.4g}; median difference '
      f'{both.rho_diff.median():.3f}; ext higher in {(both.rho_diff > 0).sum()}/{len(both)}. This is the only test.')
print('\nmixture weights (ext):')
print(e.mixture.to_string())
