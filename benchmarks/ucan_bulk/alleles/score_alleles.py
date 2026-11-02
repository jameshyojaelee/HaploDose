"""Score RNA arm allelic imbalance (arm_baf.py) against U-CAN WGS FACETS; WGD classifier; IGF2 TD comparison.

1. Arm level. Arms with >= MIN_SITES sites. Normal RNA samples set the imbalance call threshold (99th
   percentile of normal arm LLR). Per tumour: Spearman across arms of b_hat vs FACETS expected minor-allele
   fraction (purity-mixed). Per arm across patients: ROC AUC of -b_hat for FACETS LOH (LCN = 0 over >= half
   the arm) and for allelic imbalance (|TCN - 2 LCN| >= 1 over >= half the arm).
2. WGD. DNA label: patient_analysis ploidy_WGD_proxy (FACETS ploidy >= 2.7, as on CRC slide 14), and a
   second label, fraction of the autosomes with major CN >= 2 above 0.5. RNA features per tumour: allele
   (fraction of arms imbalanced, fraction with b_hat < 0.2, median b_hat) and dosage (bulk Haplodose arm log
   ratio sd, fraction of arms called gain, called loss), plus imbalanced arms not called loss. Logistic
   regression on standardised features, 5-fold stratified CV x 10 repeats (seed 0); AUC of out-of-fold
   scores; allele+dosage vs dosage-only by paired bootstrap over patients (2,000 draws, seed 1).
3. IGF2 TD (78) vs non-TD. RNA broad gains (bulk arm gain call on every arm of chr 2, 3, 12, 13, 20) and
   RNA WGD score, against the same from DNA. Fisher exact (prevalence) or Mann-Whitney (scores); Holm
   within each source (RNA calls 6, DNA 6, RNA continuous arm log ratio 5, Mann-Whitney). Unit = patient.
Usage: score_alleles.py ARMS_BAF_TSV BULK_ARMS_TSV OUTDIR
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ucan_spatial'))
from dna_truth import CEN, facets_arms  # noqa: E402

PATIENTS = Path('/gpfs/commons/home/jameslee/CRC/results/igf2_atlas_all/genomic/patient_analysis.tsv')
MIN_SITES = 50
BROAD = {'chr2': ['chr2p', 'chr2q'], 'chr3': ['chr3p', 'chr3q'], 'chr12': ['chr12p', 'chr12q'],
         'chr13': ['chr13q'], 'chr20': ['chr20p', 'chr20q']}
ALLELE = ['frac_imbalanced', 'frac_deep', 'median_b']
DOSAGE = ['dosage_sd', 'frac_gain', 'frac_loss']


def mcn_fraction(seg):
    s = seg[seg.chromosome.isin([f'chr{i}' for i in range(1, 23)]) & seg.lcn.notna() & seg.tcn.notna()]
    w = s.end - s.start + 1
    return float(np.average((s.tcn - s.lcn).ge(2), weights=w))


def cv_scores(X, y, seed=0):
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
    score, count = np.zeros(len(y)), np.zeros(len(y))
    for tr, te in RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=seed).split(X, y):
        score[te] += model.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        count[te] += 1
    return score / count


def main(arms_baf, bulk_arms, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'pip_freeze.txt').write_text(subprocess.run([sys.executable, '-m', 'pip', 'freeze'],
                                                       capture_output=True, text=True).stdout)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    baf = pd.read_csv(arms_baf, sep='\t')
    baf = baf[baf.n_sites >= MIN_SITES]
    normal = baf[baf['sample'].str.endswith('.N')]
    thr = float(normal.llr.quantile(0.99))
    baf['imbalanced'] = baf.llr > thr
    tum = baf[baf['sample'].str.endswith('.T')]

    rows, pat, skipped = [], [], []
    for s in sorted(tum['sample'].unique()):
        p = s.split('.')[2]
        try:
            dna, seg, purity, ploidy, _ = facets_arms(p, cen)
        except (FileNotFoundError, ValueError) as e:
            skipped.append(dict(sample=s, reason=str(e)[:80]))
            continue
        rows.append(dna.reset_index().assign(sample=s))
        pat.append(dict(sample=s, patient_id=p, purity=purity, facets_ploidy=ploidy, mcn2_fraction=mcn_fraction(seg)))
    dna = pd.concat(rows)
    pat = pd.DataFrame(pat)
    pd.DataFrame(skipped, columns=['sample', 'reason']).to_csv(out / 'skipped.tsv', sep='\t', index=False)
    t = tum.merge(dna, on=['sample', 'arm'], how='inner')
    t['dna_loh'] = t.loh_fraction >= 0.5
    t['dna_imbalanced'] = t.imbalance_fraction >= 0.5
    t.to_csv(out / 'arm_table.tsv.gz', sep='\t', index=False)

    # 1. Arm level.
    per = t.groupby('sample').apply(lambda g: spearmanr(g.b_hat, g.minor_expected, nan_policy='omit')[0],
                                    include_groups=False).rename('rho').reset_index()
    per.to_csv(out / 'per_sample_rho.tsv', sep='\t', index=False)
    arm_rows = []
    for arm, g in t.groupby('arm'):
        r = dict(arm=arm, n=len(g), n_loh=int(g.dna_loh.sum()), n_imbalanced=int(g.dna_imbalanced.sum()),
                 rho=spearmanr(g.b_hat, g.minor_expected, nan_policy='omit')[0])
        for col in ('dna_loh', 'dna_imbalanced'):
            y = g[col]
            r[f'auc_{col}'] = roc_auc_score(y, -g.b_hat) if 10 <= y.sum() <= len(y) - 10 else np.nan
        arm_rows.append(r)
    arms = pd.DataFrame(arm_rows)
    arms.to_csv(out / 'per_arm.tsv', sep='\t', index=False)
    calls = pd.crosstab(t.imbalanced, t.dna_imbalanced).rename(index=str, columns=str)

    # 2. WGD.
    bulk = pd.read_csv(bulk_arms, sep='\t')
    bulk = bulk[bulk['sample'].str.endswith('.T')]
    dose = bulk.groupby('sample').agg(dosage_sd=('median_log_ratio', 'std'),
                                      frac_gain=('call', lambda c: c.eq('gain').mean()),
                                      frac_loss=('call', lambda c: c.eq('loss').mean()))
    m = tum.merge(bulk[['sample', 'arm', 'call']], on=['sample', 'arm'], how='left')
    allele = m.groupby('sample').agg(frac_imbalanced=('imbalanced', 'mean'),
                                     frac_deep=('b_hat', lambda b: b.lt(0.2).mean()),
                                     median_b=('b_hat', 'median'))
    allele['frac_imbalanced_not_loss'] = m.assign(x=m.imbalanced & m.call.ne('loss')).groupby('sample').x.mean()
    lab = pd.read_csv(PATIENTS, sep='\t', usecols=['patient_id', 'ploidy_WGD_proxy', 'igf2_td'])
    d = pat.merge(lab, on='patient_id', how='left').set_index('sample').join(allele).join(dose).dropna(
        subset=ALLELE + DOSAGE + ['ploidy_WGD_proxy'])
    d['wgd_mcn'] = (d.mcn2_fraction > 0.5).astype(int)
    sets = {'dosage': DOSAGE, 'allele': ALLELE + ['frac_imbalanced_not_loss'],
            'both': DOSAGE + ALLELE + ['frac_imbalanced_not_loss']}
    wgd = {}
    for label in ('ploidy_WGD_proxy', 'wgd_mcn'):
        y = d[label].astype(int).to_numpy()
        res = {'n': int(len(y)), 'n_wgd': int(y.sum())}
        for name, cols in sets.items():
            d[f'score_{name}_{label}'] = cv_scores(d[cols].to_numpy(), y)
            res[f'auc_{name}'] = roc_auc_score(y, d[f'score_{name}_{label}'])
        for col in ALLELE + DOSAGE + ['frac_imbalanced_not_loss']:
            res[f'auc_single_{col}'] = roc_auc_score(y, d[col])
        rng = np.random.default_rng(1)
        diff = []
        for _ in range(2000):
            i = rng.integers(0, len(y), len(y))
            if y[i].min() == y[i].max():
                continue
            diff.append(roc_auc_score(y[i], d[f'score_both_{label}'].to_numpy()[i])
                        - roc_auc_score(y[i], d[f'score_dosage_{label}'].to_numpy()[i]))
        res['auc_both_minus_dosage_ci95'] = [float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))]
        wgd[label] = res
    d.to_csv(out / 'patient_table.tsv', sep='\t')

    # 3. IGF2 TD vs non-TD.
    gain_rna = bulk.assign(g=bulk.call.eq('gain')).pivot(index='sample', columns='arm', values='g')
    gain_dna = dna.assign(g=dna.dna_state.eq('gain')).pivot(index='sample', columns='arm', values='g')
    e = d[d.igf2_td.notna()]
    td = e.igf2_td.eq(1)
    tests = []
    for source, gain, wgd_score in (('rna', gain_rna, e['score_both_ploidy_WGD_proxy']), ('dna', gain_dna, e.facets_ploidy)):
        for chrom, arm_list in BROAD.items():
            x = gain.reindex(e.index)[arm_list].fillna(False).all(axis=1)
            tab = [[int((x & td).sum()), int((~x & td).sum())], [int((x & ~td).sum()), int((~x & ~td).sum())]]
            odds, p = fisher_exact(tab)
            tests.append(dict(source=source, test=f'broad_gain_{chrom}', td_rate=x[td].mean(), non_td_rate=x[~td].mean(),
                              effect=odds, p=p))
        p = mannwhitneyu(wgd_score[td], wgd_score[~td]).pvalue
        tests.append(dict(source=source, test='wgd_score' if source == 'rna' else 'facets_ploidy',
                          td_rate=wgd_score[td].median(), non_td_rate=wgd_score[~td].median(), effect=np.nan, p=p))
    # Arm gain calls (|z| > 3) are conservative, so RNA is also tested on the continuous arm log ratio (mean over
    # the chromosome's arms; medians reported in the rate columns), Mann-Whitney, its own Holm family of 5.
    lr = bulk.pivot(index='sample', columns='arm', values='median_log_ratio').reindex(e.index)
    for chrom, arm_list in BROAD.items():
        x = lr[arm_list].mean(axis=1)
        tests.append(dict(source='rna_log_ratio', test=f'arm_log_ratio_{chrom}', td_rate=x[td].median(),
                          non_td_rate=x[~td].median(), effect=np.nan, p=mannwhitneyu(x[td], x[~td]).pvalue))
    tests = pd.DataFrame(tests)
    tests['p_holm'] = tests.groupby('source').p.transform(lambda p: multipletests(p, method='holm')[1])
    tests.to_csv(out / 'td_tests.tsv', sep='\t', index=False)

    summary = dict(n_tumours=int(t['sample'].nunique()), n_skipped=len(skipped), normal_llr_threshold=thr,
                   n_normal_arms=int(len(normal)),
                   median_per_tumour_rho=float(per.rho.median()),
                   median_arm_auc_loh=float(arms.auc_dna_loh.median()),
                   median_arm_auc_imbalance=float(arms.auc_dna_imbalanced.median()),
                   imbalance_call_vs_dna=calls.to_dict(), wgd=wgd, n_td=int(td.sum()), n_non_td=int((~td).sum()),
                   correction='Holm within source (rna calls 6, dna 6, rna_log_ratio 5) for TD tests; unit = patient')
    (out / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))
    print(tests.round(4).to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
