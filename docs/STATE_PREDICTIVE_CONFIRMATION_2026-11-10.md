# A predictive state safeguard fails recovery, so reject it

The implemented discovery-to-confirmation predictive comparison does not solve
the clone false-split/power tradeoff. It retains original case20261107 but rejects
original case20261109; three of four known false splits favor continuous state,
and the fourth is numerically unresolved. The prespecified six-case screen
fails. Reject expansion and do not tune its threshold, override the original
recovery assertion or claim5% error control. Production and frozen results are
unchanged. All eight CPU jobs complete0:0.

## Correction tested, rather than another computation-only review

The preceding discrete-versus-continuous BIC comparison used discovery counts.
This experiment fits those two distributions on discovery only, then compares
their predictions of confirmation counts conditional on the SAME cell's latent
regional state. Its target is an added safeguard against transcription-induced
subdivisions while preserving all original recoveries. Neither distribution is
assumed to be a certified DNA genotype model.

[Executable implementation](../benchmarks/hardening/state_predictive_confirmation.py),
[design frozen before outcomes](../dev/state_predictive_confirmation_20261110/DESIGN.md),
[exact case selection](../dev/state_predictive_confirmation_20261110/CASES.json).
Select both difficult original cases whose preceding BIC safeguard failed, plus
ALL four false splits in the earlier phase6 candidate pilot. Each has exactly
one accepted chr12 event. Keep its cells, genes, bins, discovery offsets,
reference-only NB sizes, support mask and proposal. Do not choose a more favorable
region, parent, sample or count allocation after seeing predictive scores.

Let D and C be discovery and confirmation counts. The statistic is
log p(C|D,mixture)-log p(C|D,continuous), evaluated as joint(D,C) minus discovery(D)
integrals. Both integrals share one latent state per cell; separately integrating
two independent states would discard membership information. The alternative
has two Gaussian state components with shared variance; the comparator has one.
Use the corrected variance optimizer and exact zero-variance boundary, with
unchanged starts/constraints. All24/48/96-order model fits freeze before any
confirmation-count array is decoded. Reference expectations and upstream
selection are conditional, not established independent by this new fit.

Observation assumption: conditional independent NB genes and an adequate
shared-state distribution. Beta-binomial count folds are independent conditional
on the state when the NB size law is correct; estimated or misspecified sizes
need not preserve that property. This experiment does not certify that assumption.
Regional transcription/reference variation can have the same conditional RNA
law as a CNA population. Independent corresponding DNA, ideally same-cell for
membership, is the observation that can distinguish genomic interpretations.

No new scientific cutoff is fitted. Positive predictive advantage is a necessary
bounded-screen condition for retaining each true recovery; negative advantage
is required for each known false split. Both must be resolved beyond the fixed
quadrature/rounding diagnostic. Passing would permit subsequent full original
recovery and fresh complete-workflow calibration; it would not itself establish
an error rate. Failure rejects expansion. The one-sided95% upper specimen
false-split bound<=5%, simultaneous primary assessment and power requirements
remain unchanged.

## Independent numerical result

An exact Fraction calculation of integer NB probabilities and component weights
gives predictive log score-0.43183625835986644. Error2.78e-16 is below the
operation-derived2.18e-11 budget. Independently integrate actual joint NB PMFs
with SciPy adaptive quadrature, including constants, then divide by the
discovery density. Error5.00e-16 is below integration-error plus rounding budget
5.18e-11. Integration beyond12 standard deviations is bounded explicitly by
twice the normal survival probability because each PMF is<=1. Component swap
also passes. [Numerical report](../dev/state_predictive_confirmation_20261110/numerics_v1/report.json).

These references establish computational agreement for the tested integrals,
not inferential calibration. The real-case two-times successive-order diagnostic
is explicitly a convergence heuristic, not a proven integration bound or an
interval with biological coverage. An unresolved sign stays unresolved.

## Complete fixed-case result

| Case | Cells | Genes | Conditional mixture log-score advantage | Numerical diagnostic | Resolved | Proposed guard |
|---|---:|---:|---:|---:|---|---|
| Original20261107 |512|93|4.561159872|0.000850044|yes|retain|
| Original20261109 |502|93|-3.923530592|0.000526640|yes|reject, fails recovery|
| State202611096001 |599|92|-0.000012977|0.001622507|no|abstain|
| State202611096002 |600|94|-2.925425639|0.000704574|yes|reject|
| State202611096004 |599|93|-5.748000004|0.001315536|yes|reject|
| State202611096006 |600|94|-10.047175953|0.001285647|yes|reject|

These are selected observations, not600-cell biological replications. Biological
n=0. All four false-split cases were already exposed development data, so their
rejection is not a prospective false-split estimate. Missing selected cells are
not silently promoted to a complete-specimen assignment accuracy denominator.

[All order-specific scores and model fits](../dev/state_predictive_confirmation_20261110/summary_v1/report.json).
Each worker saves discovery parameters before confirmation, exact cell/gene
identities and per-cell predictive scores/memberships. Full counts equal
discovery plus confirmation exactly. No true clone labels enter fitting.

[Observed pattern, not a guarantee] A confirmation model that retains membership
information still favors continuous state on one true simulated recovery case.
This falsifies this safeguard as a repair satisfying the original cases. It
does not prove that the case is universally unrecoverable, that RNA contains no
clone information, or that any real tumour has a continuous rather than genomic
origin. A changed threshold after these outcomes would not resolve that issue.

## Compute, changes and preservation

All jobs use CPU/v6&amd, the read-only wheel validated by431 Python/4 Rust tests,
neutral directory, protected output paths and recorded import/artifact hashes.
No arrays or shared installs. [Exact commands](../dev/state_predictive_confirmation_20261110/jobs.json).

| Job | Role | Account/QOS | CPUs/RAM/time | Outcome |
|---|---|---|---|---|
|23155390|Independent numerical references|nslab|1/4GB/1h|complete|
|23155391|Original20261107|innovation|1/8GB/2h|complete|
|23155392|Original20261109|nslab|1/8GB/2h|complete|
|23155393|State202611096001|innovation|1/8GB/2h|complete|
|23155394|State202611096002|nslab|1/8GB/2h|complete|
|23155395|State202611096004|innovation|1/8GB/2h|complete|
|23155396|State202611096006|nslab|1/8GB/2h|complete|
|23155397|Gated summary|nslab|1/4GB/1h|complete|

The only analytical source addition is the research implementation; no production
state/clone code, original test, scientific threshold, previous claim entry or
frozen result changes. The active goal remains incomplete. The next practical
profile step makes the preceding successful RNA stack usable from identified
count matrices while retaining its proven estimator and explicit limitations.
