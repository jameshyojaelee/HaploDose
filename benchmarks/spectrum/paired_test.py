#!/usr/bin/env python3
"""SPECTRUM claim 4 primary test (docs/STATUS.md, claim 4), over patients.

Input: PATIENTS (TSV with a 'patient' column: the prespecified eligible patients, fixed from the DNA clone tables
before any RNA tool ran) and SCORES_ROOT/<patient>/<tool>/scores.tsv written by score.py (level=patient row).
Tools: haplodose and the REQUIRED comparators below. Amendment 2026-11-03 (before any CNV tool ran on SPECTRUM
scRNA; docs/STATUS.md "Claim 4 SPECTRUM run rules"): the family is fixed at the 4 required comparators whether or
not they produce output; the earlier rule (family = comparators that ran) let a missing comparator drop out.
XClone and CalicoST are not installed and are not part of the claim. Per patient, a missing or no_output score
counts as AUROC 0.5 (no information); completion is reported per tool alongside the test. Because a comparator
that fails everywhere would sit at 0.5 and be beaten trivially, SOTA also requires every required comparator to
have scored output for at least half the patients (MIN_COMPLETION); otherwise the claim is 'not established'.
Descriptive sensitivity: the same tests over patients where both tools scored (no 0.5 fill).
Descriptive sensitivity (review F5, added before any output): common arms. Per patient and pair, AUROC over the
arms both tools measured (arm_scores.tsv 'measured', patient row), over patients where it is defined for both.

Test, HaploDose vs each comparator in the family: paired two-sided Wilcoxon signed-rank over patients
(benchmarks/heldout/analyse.py wilcoxon_p: zero differences dropped; p NaN if every difference is 0), Holm over
the family (analyse.py family: an untestable comparison enters with p = 1). Mean paired difference with a 95%
percentile bootstrap CI (10,000 resamples of patients, seed 20261103, a fresh generator per comparison so a
comparison does not depend on which others ran). Superior: Holm p < 0.05 and CI lower bound > 0.
SOTA (claim 4) = superior to every comparator in the family.

Usage: paired_test.py PATIENTS SCORES_ROOT OUT
Writes OUT/per_patient.tsv and OUT/primary.json (OUT must not exist or be empty).
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

H = Path(__file__).resolve().parents[2]
SEED, B = 20261103, 10_000
TOOL = 'haplodose'
COMPARATORS = ['infercnvpy', 'copykat', 'scevan', 'numbat']   # required family, fixed before any run
MIN_COMPLETION = 0.5

_spec = importlib.util.spec_from_file_location('heldout_analyse', H / 'benchmarks/heldout/analyse.py')
analyse = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(analyse)


def common_arm_auroc(root, patient, a, b):
    """(AUROC a, AUROC b) over arms both tools measured; NaNs when undefined or no arm table."""
    t = {}
    for tool in (a, b):
        f = Path(root) / patient / tool / 'arm_scores.tsv'
        if not f.exists():
            return np.nan, np.nan
        x = pd.read_csv(f, sep='\t')
        x = x[x.specimen.eq('patient')].set_index('arm')
        if 'measured' not in x:
            return np.nan, np.nan
        t[tool] = x
    arms = t[a].index[t[a].measured.astype(bool)].intersection(t[b].index[t[b].measured.astype(bool)])
    y = t[a].truth.reindex(arms)
    if not 0 < y.sum() < len(y):
        return np.nan, np.nan
    return tuple(float(roc_auc_score(y, t[k].score.reindex(arms))) for k in (a, b))


def patient_auroc(root, patient, tool):
    """(AUROC or NaN, status) from score.py's patient row."""
    f = Path(root) / patient / tool / 'scores.tsv'
    if not f.exists():
        return np.nan, 'missing'
    t = pd.read_csv(f, sep='\t')
    r = t[t.level.eq('patient')]
    if len(r) != 1 or r.status.iloc[0] != 'ok' or pd.isna(r.diff_auroc.iloc[0]):
        return np.nan, str(r.status.iloc[0]) if len(r) == 1 else 'unreadable'
    return float(r.diff_auroc.iloc[0]), 'ok'


def compare(a, b):
    d = (a - b).to_numpy(float)
    rng = np.random.default_rng(SEED)
    boot = d[rng.integers(0, len(d), (B, len(d)))].mean(axis=1)
    ci = [float(x) for x in np.quantile(boot, [0.025, 0.975])]
    return dict(n=len(d), mean_diff=float(d.mean()), median_diff=float(np.median(d)), ci95=ci,
                a_higher=int((d > 0).sum()), b_higher=int((d < 0).sum()), p=analyse.wilcoxon_p(d))


def main(patients, root, out):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    pts = pd.read_csv(patients, sep='\t').patient.astype(str).tolist()
    rows = []
    for p in pts:
        for tool in [TOOL] + COMPARATORS:
            v, s = patient_auroc(root, p, tool)
            rows.append(dict(patient=p, tool=tool, auroc_raw=v, status=s))
    per = pd.DataFrame(rows)
    per['auroc'] = per.auroc_raw.fillna(0.5)
    per.to_csv(out / 'per_patient.tsv', sep='\t', index=False)
    wide = per.pivot(index='patient', columns='tool', values='auroc').reindex(pts)
    ok = per[per.status.eq('ok')].groupby('tool').size().reindex([TOOL] + COMPARATORS, fill_value=0)
    tests = analyse.family([dict(comparator=c, **compare(wide[TOOL], wide[c])) for c in COMPARATORS])
    raw = per.pivot(index='patient', columns='tool', values='auroc_raw').reindex(pts)
    both = []
    for c in COMPARATORS:
        m = raw[TOOL].notna() & raw[c].notna()
        both.append(dict(comparator=c, **(compare(raw[TOOL][m], raw[c][m]) if m.sum() else dict(n=0, p=np.nan, ci95=[np.nan, np.nan]))))
    both = analyse.family(both)
    complete = {c: bool(ok[c] >= MIN_COMPLETION * len(pts)) for c in COMPARATORS}
    common = []
    for c in COMPARATORS:
        v = pd.DataFrame([common_arm_auroc(root, p, TOOL, c) for p in pts], columns=['a', 'b']).dropna()
        common.append(dict(comparator=c, **(compare(v.a, v.b) if len(v) else dict(n=0, p=np.nan, ci95=[np.nan, np.nan]))))
    common = analyse.family(common)
    summary = dict(
        n_patients=len(pts), unit='patient', endpoint='clone-difference AUROC (score.py level=patient)',
        family=COMPARATORS, no_output_any_patient=[c for c in COMPARATORS if ok[c] == 0],
        scored_patients=ok.to_dict(), failures_scored_0_5=(len(pts) - ok).to_dict(),
        mean_auroc={t: float(wide[t].mean()) for t in [TOOL] + COMPARATORS},
        tests=tests, completion_ok=complete,
        sota=all(t['superior'] for t in tests) and all(complete.values()),
        sensitivity_both_scored=both, sensitivity_common_arms=common,
        method=('paired two-sided Wilcoxon over patients, Holm over the 4 required comparators; bootstrap 95% CI of '
                f'the mean difference ({B} patient resamples, seed {SEED}); missing patient AUROC = 0.5'))
    (out / 'primary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main(*sys.argv[1:])
