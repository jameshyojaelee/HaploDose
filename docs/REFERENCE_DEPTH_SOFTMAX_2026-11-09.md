# Continuous depth regression: selective benefit, no general default

Continuous reference-cell regression with heldout-selected shrinkage reduces
noise under zero slopes and improves the small-reference binned simulation.
It does not solve narrow depth support or a reference-specific programme. Its
benefit does not transfer to a gene-level fixture, where simple pooling has
slightly lower RNA-profile MSE. Retain the prototype and negative evidence;
do not promote a general correction or claim improved clone/DNA performance.

## Estimand, model and identifying assumption

The first experiment targets expected normalized pooled RNA abundance at the
query exposure distribution. It reuses all48 previously exposed scenarios from
the [estimated-slope panel](DEPTH_REFERENCE_ESTIMATION_2026-11-09.md), with200 fit
references except30 in the small condition,100 heldout references and100 query
exposures. Eight simulation replicates per condition, biological n=0. Query
exposures are known simulation inputs, not the observed library totals available
to a deployed inference workflow.

The mean model is softmax(b+gamma*logL). Fit intercepts and slopes to individual
reference counts with count-weighted cross entropy, optionally penalizing the
squared centred slopes. This is a composite mean-model objective for the
gamma-Poisson generated counts; their conditional joint law is not asserted
multinomial. Its Hessian is not a biological uncertainty estimate. Observed
library totals also differ from the generating exposures.

Reference depth association must transport to the query. A reference-specific
expression programme is the strongest competing explanation, represented
explicitly in the confounded challenge. Better prediction of heldout references
cannot distinguish a transferable nuisance effect from that programme. Real
interpretation needs independent reference metadata, depth-matched RNA
replication and corresponding DNA at the claimed resolution.

## Frozen comparison

[Design](../dev/reference_depth_softmax_20261109/DESIGN.md) and
[implementation](../benchmarks/hardening/reference_depth_softmax.py) define eleven
candidates: pooled zero slopes, unpenalized regression and nine ridge strengths
scaled to the reference likelihood curvature. Scaling makes penalties interpretable
relative to the data objective, not biological thresholds. The last intercept
and slope are anchored; the penalty uses all centred slopes. Standardizing
log-depth is purely numerical.

Fit using reference fit counts only; select minimum heldout Poisson deviance,
with a fixed simplicity tie rule. The holdout is development selection data,
not clone confirmation. Predictions, scores and choices are hashed before
generating targets are opened. Every original mechanism condition and comparator
is retained. Neither biological thresholds nor scientific criteria changed.

Mean null centred log-rate MSE across the eight replicates:

| Condition | Current expectation | Quintile/intercept correction | Selected continuous regression | Simple pooled reference |
|---|---:|---:|---:|---:|
| Supported |0.00399348|0.00026261|0.00026378|0.00075295|
| Extrapolated |0.30930058|0.00639593|0.00704459|0.32397516|
| Narrow |0.38793159|0.62396112|0.38272043|0.39259781|
| Small |0.33021826|0.33021826|0.04043890|0.33021826|
| Zero slopes |0.00033397|0.00580445|0.00033634|0.00030965|
| Reference-specific programme |0.33685558|0.18840044|0.19130314|0.34521562|

Selection chooses pooled rates in5/8 zero-slope and4/8 narrow cases; the remaining
cases use strong shrinkage. Small-reference cases use weak shrinkage or no
penalty, avoiding the quintile estimator's60-cell support cutoff. This conditional
benefit does not justify reducing a scientific support requirement.

The selected mean first-region gain remains negative under narrow support
(−0.62466) and confounding (−0.61655), despite a generating positive RNA effect
0.24328. A general correction/gain-preservation claim therefore fails. Small
reference mean gain is0.10819, not complete signal recovery. All individual
outcomes, candidate scores and chosen complexities are saved.

## Independent numerical and leakage evidence

[Audit](../benchmarks/hardening/reference_depth_softmax_audit.py) checks21120
predicted bins,528 scalar heldout scores and1440 target metric values. It
independently checks score equations, positive local Hessians, six directional
Hessian calculations at two complex steps, and bin permutation/anchor invariance.
Maximum directional discrepancy is2.78e-17; anchor permutation log-prediction
change is1.11e-15. Operation-count precision bounds govern numerical comparisons,
not biological error rates.

Coarse/fine optimization controls produce no case where the winning heldout score
gap overlaps the observed score-refinement discrepancy. This is a numerical
resolution check, not statistical certainty about a model ranking. Solver flags,
stationarity residuals and curvature remain recorded separately; a stationary
working-model optimum is not scientific validation.

A complete fitting replay with altered generating targets preserves every
coefficient, predicted mean, heldout score and selection exactly, while evaluated
errors change. Generating truth does not select the model. The mutation is an
information-leakage test, not independent biological evidence.

## Gene-level aggregation challenge

[Separate prespecified design](../dev/reference_depth_softmax_20261109/GENE_FIXTURE_DESIGN.md)
uses the existing raw-count parity fixture, known generating type0 identity,
30 fit and30 heldout reference cells,160 query cells and48 bins. Counts select
and freeze the model before generating means are accessed. No query clone labels
define the target or fitted reference.

This experiment's estimand is the pooled query RNA contrast to the unaltered
type0 reference counterfactual at the query's true exposure distribution. It
differs from the earlier unstandardized pooled-RNA target; scores must not be
compared across those targets as a method improvement. Reference identity and
counterfactual depth exposures are oracle simulation information, not established
from biological RNA alone.

For truth construction only, two noiseless reference cells identify relative
gene slopes/intercepts. Reconstructing all60 reference gene-mean vectors agrees
to2.66e-15 maximum log error. The law then generates the unaltered reference
counterfactual at query exposures before exact gene-to-bin aggregation. Count
fitting still uses observed library totals. Heterogeneous gene slopes need not
aggregate into the fitted bin-level softmax law, providing a distinct challenge.

| Method | Unsmooth MSE | Unsmooth Pearson | Unsmooth Spearman | Shared-partition MSE |
|---|---:|---:|---:|---:|
| Pooled reference |0.00341151|0.923648|0.905667|0.00756472|
| Quintile/intercept |0.00341151|0.923648|0.905667|0.00756472|
| Unpenalized continuous |0.00345025|0.922964|0.905558|0.00758071|
| Heldout-selected continuous |0.00341587|0.923661|0.905015|0.00756798|

The selected ridge strength is the curvature-scale candidate. It has slightly
higher Pearson but higher MSE and lower Spearman than pooling; no dominance is
claimed. The identical previously frozen RNA-selected partition worsens MSE for
every method at this RNA target. It is not DNA segmentation truth. These are
descriptive results from one exposed fixture, without biological uncertainty.

[Gene audit](../benchmarks/hardening/reference_depth_gene_audit.py) independently
reconstructs14400 reference gene means,38400 query counterfactual means,48 target
bins and24 metric values using scalar pooling, ranks and correlations. All pass.
No universal transfer from the simpler binned simulation is demonstrated.

## Compute, preservation and decision

| Job | Purpose | Account | Exit | Runtime |
|---|---|---|---:|---|
|23083219|48-case candidate family|innovation|0|27s|
|23083646|Independent numerical audit and truth replay|nslab|0|39s|
|23084095|Gene-level fixture|innovation|0|6s|
|23084253|Independent gene/metric audit|nslab|0|17s|

All jobs used allocated CPU nodes, immutable source snapshots, neutral working
directories and the validated noneditable wheel. Resources and identities are
in [jobs.json](../dev/reference_depth_softmax_20261109/jobs.json). Existing package,
native, test and README files are unchanged; historical artifacts remain intact.

Retain simple pooling as the baseline and do not add this general default.
No large bootstrap expansion is justified by the failed transfer/generalization
checks. Any future mode-specific use needs reference uncertainty including
complexity selection, support/abstention behavior, consistent reference-depth
transport through the entire inference model and frozen complete-workflow
false-split/power tests. Unidentified reference programmes need independent
information, not a more elaborate optimizer. Original three recovery failures,
whole-specimen clone error/power and biological validation remain unresolved.
