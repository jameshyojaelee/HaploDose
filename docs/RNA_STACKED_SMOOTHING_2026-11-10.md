# An implemented RNA profile improvement, with a failed bulk transfer

The fixed conditional RNA stack improves 1 Mb DNA-profile Pearson from
**0.667955 to 0.707226** against width-five pseudobulk across all181 previously
exposed Kinker lines, with177 wins and4 losses. It also beats the preceding
predictive partition mean0.685905. All five prespecified conditional comparison
tests pass their simultaneous interval and Holm criteria. This is retrospective
regression evidence, not independent biological validation or universal
superiority. The same stack fails the bulk replacement gate and is rejected
for that mode. Frozen results, old comparator scores and scientific criteria
remain unchanged.

## What was implemented

[Estimator and runners](../benchmarks/hardening/rna_stacked_smoothing.py),
[bulk runner](../benchmarks/hardening/bulk_stacked_smoothing.py),
[fixed comparison](../benchmarks/hardening/rna_stacked_summary.py),
[initial design](../dev/rna_stacked_smoothing_20261110/DESIGN.md),
[181-line/bulk design](../dev/rna_stacked_smoothing_20261110/REGRESSION_DESIGN.md).

Estimand: a regional descriptive natural-log RNA abundance contrast relative to
the declared reference. Inputs are the existing prepared gene contrasts, gene
order, bins and contiguous arm/gap blocks. Kinker counts were reconstructed from
CPM, so this experiment does not claim a raw-count sampling likelihood. There
is no Poisson thinning, transferred bulk count model or new DNA calibration.

The expert bank contains the existing geometric running-median widths, plus
the existing piecewise-mean Gaussian working smoother. That smoother retains
its robust adjacent-difference scale and3log(block length) penalty. No new
regularization constant, DNA-selected width, reference or temperature is added.
Fold assignment is the existing within-bin gene rank modulo ten. Each eligible
target is removed from its base predictor's training data; singleton genes
remain training observations and are never prediction targets. Every eligible
gene receives exactly one out-of-training base prediction.

Fit nonnegative weights summing to one by minimizing absolute prepared-RNA
prediction error. The weights use those targets: the resulting objective is
**not an independent confirmation score or unbiased cross-validation error**.
Upstream gene selection, reference fitting and normalization are conditioned
on. Full-data expert predictions are combined using the fitted weights.

The sparse SciPy/HiGHS solution is checked independently using the feasible
clipped L1 dual. A conservative rounding budget scales with operation count,
machine epsilon and data magnitude. The primal-minus-dual gap bounds numerical
optimization uncertainty; it is not a statistical interval. Use the stack only
if its risk plus the gap and rounding budget is below every single expert.
Otherwise retain the installed median selector; all-singleton input explicitly
abstains. A nonidentified positive segmentation scale excludes that expert.

Identifying assumption for a DNA interpretation: the reference contrast and
regional expression response reflect genomic dosage at the scored resolution.
The strongest alternative is regional transcriptional regulation/reference
composition producing the same RNA shape. Smoothing cannot distinguish these.
Corresponding independently processed DNA, replication across biological
specimens and prospective locked comparisons are needed for stronger claims.

## Independent numerical evidence and locked cases

[Numerical report](../dev/rna_stacked_smoothing_20261110/numerics_v2/report.json)
checks an exact Fraction-based two-expert L1 optimum1/3 at weights1/2,1/2,
independently enumerated breakpoints, feasible dual, row permutation, target
perturbation invariance of base predictions, constant fallback and singleton
abstention. No universal epsilon is used to decide scientific success.

Before broader scoring, both locked diagnostic cases improve against all five
comparators: `2313287_STOMACH` and `42MGBA_CENTRAL_NERVOUS_SYSTEM`.
Their meanPearson is0.819845 versus retained predictive0.793928 and sm5
0.762426. The RNA objective improves over the best single expert by
0.003926/0.002339; feasible-dual gaps are below6e-9. All10 previous scores
replay exactly. [Locked report](../dev/rna_stacked_smoothing_20261110/locked_v2/report.json).

Source_v1 was not executed. Source_v2 fixes HiGHS to one thread through the
installed SciPy-supported pass-through option; its visible OptimizeWarning is
retained. Source_v3 adds chunked runners only; the estimator functions are
unchanged fromv2. No unsuccessful scientific candidate was hidden by this
versioning.

## All181 regression lines, no case removal

The original order, genes, reference membership/UMI DNA weights, windows and
relative DepMap contrast are retained. All905 previous scores are exact.
Predictions and RNA choices freeze before DNA scoring. All181 lines have
defined results; shared query/reference identities yield eight connected
resampling groups, not181 independent reference panels.

| Method | Mean 1 Mb Pearson |
|---|---:|
| Raw bin pseudobulk |0.5959399404|
| Width-five pseudobulk |0.6679553826|
| Pseudobulk MDL |0.6683618390|
| Retained predictive median |0.6802181879|
| Preceding predictive partition |0.6859045160|
| Conditional RNA stack |**0.7072257053**|

Fixed primary familyfive,20,000 whole-group draws, seed2026111001; simultaneous
basic Bonferroni intervals and Holm p<.05 are both required. These intervals
are approximate conditional bootstrap summaries with only eight groups and
fixed upstream fits. Empirical interval coverage was not established here.

| Comparator | Mean paired gain | Wins/losses | Simultaneous interval | Holm p |
|---|---:|---:|---|---:|
| Raw bin pseudobulk |0.111285765|177/4|[0.096078930,0.128452420]|0.000249988|
| Width-five pseudobulk |0.039270323|177/4|[0.032247085,0.045814846]|0.000249988|
| MDL pseudobulk |0.038863866|154/27|[0.031779663,0.046370828]|0.000249988|
| Retained predictive median |0.027007517|160/21|[0.021588939,0.032113914]|0.000249988|
| Preceding predictive partition |0.021321189|151/30|[0.015668687,0.027192782]|0.000249988|

[Complete scores and fixed-family report](../dev/rna_stacked_smoothing_20261110/summary_v3/report.json),
[all paired profiles](../dev/rna_stacked_smoothing_20261110/summary_v3/paired_profiles.tsv).
This establishes a task-specific regression improvement, not a new unopened
evaluation. Relative DepMap measurements remain uncertain cross-culture
reference measurements, not same-cell clone labels.

## Bulk transfer fails, so do not replace its better estimator

Use the same fixed stack on all48 preceding raw-bulk prepared contexts, with
the same20 normal references,24/24 patient split, normal centering, masks,
loss threshold and all39 old methods. All39 old score tables replay exactly;
the candidate is the40th method. Undefined endpoints remain visible.

On the comparison side, correlation is defined for23/24 patients and loss
recall for13. MeanPearson of the stack is0.7035111344, below strongest
pseudobulk0.7302613671 and posterior-averaged reference median0.7648167082.
It wins9/23 and2/23 respectively. Its loss recall0.538415751 includes FPR
0.049186284, compared with pseudobulk recall0.399679487/FPR0.041244451 and
posterior averaging recall0.602747253/FPR0.043737174.

The prespecified positive-gain gate fails. Reject bulk replacement and skip
further bulk uncertainty for this losing candidate. Keep the better existing
bulk path and its earlier failed69-family superiority decision.
[Bulk report](../dev/rna_stacked_smoothing_20261110/bulk_v3/report.json),
[all40 methods/endpoints](../dev/rna_stacked_smoothing_20261110/bulk_v3/summary.tsv).

## Compute and reproducibility

Eight jobs complete0:0, all CPU/v6&amd, isolated validated read-only wheel,
neutral working directory. No arrays, shared installs or frozen output writes.

| Job | Role | Account/QOS | CPUs/RAM/time | Outcome |
|---|---|---|---|---|
|23154535|Independent numerics|nslab|1/4GB/1h|complete|
|23154536|Locked two cases|innovation|1/8GB/2h|complete|
|23154686|Lines0:46|innovation|1/8GB/2h|complete|
|23154687|Lines46:92|nslab|1/8GB/2h|complete|
|23154688|Lines92:138|innovation|1/8GB/2h|complete|
|23154689|Lines138:181|nslab|1/8GB/2h|complete|
|23154690|Bulk48|innovation|1/8GB/2h|complete|
|23154691|Afterok fixed-family summary|nslab|1/4GB/1h|complete|

[Exact commands and dependencies](../dev/rna_stacked_smoothing_20261110/jobs.json).
Every worker records executed script/environment/wheel/native identities,
runtime, memory, backend, source and transformed output hashes. Initial dirty
source is preserved in the task baseline; frozen evidence is unchanged.

## Package delivery and next constraint

The proven estimator is being delivered as the opt-in
`haplodose.stacking.stacked_rna_profile` API. Its full tests, installed181-line
parity and matched field comparison are recorded separately in
[native delivery](NATIVE_STACKED_PROFILE_2026-11-10.md). No bulk or clone default
changes follow from profile correlation.

The three original recovery failures are already solved by the research gene
candidate on their original assertions, but that candidate falsely splits4/10
continuous-state nulls. Its promotion remains rejected; production remains2/5.
See [full recovery and false-split/power tradeoff](CONDITIONED_GENE_SCORE_2026-11-09.md).
This profile improvement neither repairs nor relaxes that failure, the5%
specimen confidence requirement, uncertainty coverage, or biological truth gaps.
