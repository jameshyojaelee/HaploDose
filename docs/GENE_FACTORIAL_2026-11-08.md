# Direct-gene preprocessing: no advantage over the matched simple baseline

The fixed eight-way factorial does not establish superiority over the matched
naive profile. No production default is changed. This is development evidence
from already exposed specimens, not independent biological validation.

## Measurement and estimand

The factorial varies mapped/full library depth, pooled-rate/equal-donor-log
reference construction, and half-count/one-CPM pseudocounts. All eight variants
use the same genes, common callable support, RNA-only partition and DNA reference.
Seventeen historical methods remain, giving 25 methods and 51 paired contrasts
across three endpoints (153 hypotheses, Holm correction). The correlation
estimand is the mean within-patient difference in window correlation with DNA.
There are 24 calibration and 24 comparison patients; correlation comparisons
retain 23 comparison patients and explicitly record one missing patient.

Uncertainty resamples reference donors, calibration patients and comparison
patients. Only eight saved reference resamples are available. A million draws
reduces Monte Carlo error conditional on that construction; it does not establish
coverage under biological sampling, reference adequacy, or independent replication.

## Verified results

The formula audit independently reconstructs all 72 transforms (eight variants
and nine fits). The array audit retains 171 historical arrays and checks 10,800
patient-method-refit rows. The million-draw correlation audit independently
checks 102 probability counts per checkpoint, 368 interval endpoints per
checkpoint and 714 scalar draw endpoints. Scientific and numerical families
remain 153 and 459. Loss uncertainty is still being replayed with the corrected
exact-patient-weight scorer; its unassessed precision is represented by worst-case
probability bounds in the correlation audit.

Direct pooled-gene versus intermediate pooled-median aggregation has a mean
Pearson advantage of 0.03592316. The programme-controlled counterpart has an
advantage of 0.03161900. Both retain positive simultaneous intervals and Holm
decisions after Monte Carlo bounds; the latter's conservative adjusted probability
upper bound is 0.04204521. These are matched aggregation effects, not evidence
that programmes themselves improve accuracy.

None of the eight variants improves correlation over `naive_own_partition` under
the declared family. Pearson differences range from -0.00547717 to approximately
zero; all eight adjusted probabilities are 1. This negative result is retained.

The full-depth/log-mean/CPM variant is an algebraic naive-baseline control:
`ln(counts/L + 1e-6) = ln(CPM + 1) - ln(1e6)`, and reference subtraction cancels
the constant. An independent diagnostic over 432 profiles finds no segmentation
boundary disagreement. After removing each profile's scalar translation, maximum
residual is 5.461923e-8 natural-log units, attributable to differing float32 cast
locations. Different centering can still change thresholded loss calls; this is
not a proof of loss-decision equivalence.

## Decisions and evidence

Retain direct-gene aggregation as a supported internal correction to intermediate
median aggregation. Do not promote factorial settings as a new biological
mechanism or superior method. Complete the exact loss replay and fixed-budget
loss precision audit before drawing loss conclusions. Frozen prior results remain
unchanged, including the failed first CSV-parser audit.

Artifacts are under `dev/gene_factorial_20261108`: `panel_v1`, `audit_v2`,
`formula_v1`, `precision_v1`, and `precision_audit_v1`. Compute jobs 22998578,
22998713, 22998796, 22998942 and 22998989 produced these records. Baseline
equivalence is independently recorded by job 22999976 under
`dev/naive_equivalence_20261108/diagnostic_v1`.

The full scientific objective remains unmet. Neither these correlations nor
numerically precise resampling validate cell-level clones or absolute DNA dosage.

## Complete exact-loss refinement

All25 continuation jobs23000360–23000384 complete successfully, retaining exactly
one million fixed draws including the admitted20k prefix. Collector23000385
completes in2m39s and verifies every resampling identity, gap-free chunk coverage,
2652 independently recomputed scalar full-scorer effects, interval interpolation
and the full153-test family. Independent audit23000411 completes in16s and verifies
all51 loss probability counts,188 interval endpoints and153 exact-rational Holm
decisions. Evidence is under `dev/factor_loss_precision_20261108/{complete_v1,
final_audit_v1}`. Earlier20k analyses and every source snapshot remain unchanged.

There are zero loss-recall superiority decisions. Each of the eight variants
has point loss-recall difference0 versus `naive_own_partition` and Holm-adjusted
p=1. Family intervals remain wide: for example mapped/pooled/half is
[-0.33215,0.27616]. Thirteen comparison patients contribute loss truth and11 are
recorded missing for that endpoint. The result is absence of demonstrated
improvement under the declared analysis, not established equivalence or adequate
power. Correlation positives remain four rank and nine Pearson contrasts across
the full family; they do not repair the loss-side scientific gap.

Decision: no factorial setting becomes a scientific default, and no superiority
or population-calibration claim is added. Further gains require an evidenced
measurement/model intervention and adequate independent truth, rather than
additional Monte Carlo draws on this exposed panel.
