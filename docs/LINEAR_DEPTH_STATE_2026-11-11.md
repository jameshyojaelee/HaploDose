# A lower-dimensional depth trend for membership confirmation

The goal remains accurate RNA-derived regional genomic profiles and supported
CNA-defined populations, with independent evidence at the claimed resolution.
This experiment implements a new performance candidate for the unresolved
membership/continuous-state tradeoff. It replaces eight categorical nuisance
means with a common intercept and slope in discovery-only outside-event depth.
It does not change a support, posterior, confidence or error threshold. Execution
and complete scientific results are recorded below as they become available.

## Why this change, and what would falsify it

Earlier fixed/soft prevalence models lose original20261109. Balanced prevalence
can restore its conditional predictive sign but accept state202611096004. A
separate shared-prevalence model passes a smaller development screen but later
loses fresh recovery and false-split criteria. These are different retained
experiments, not interchangeable success claims:
[depth-conditioned model](DEPTH_STATE_PREDICTIVE_2026-11-10.md),
[prevalence tradeoff](PROFILE_MEMBERSHIP_PRIOR_2026-11-10.md),
[shared-prevalence mechanism](POOLED_PREVALENCE_2026-11-10.md),
[complete-workflow pilot](POOLED_PREVALENCE_PILOT_2026-11-10.md).

[Inference] Freely fitted locations in eight strata may lose information or
exaggerate conditional separation. That hypothesis is not established by the
earlier failures. The smallest new analysis tests one lower-dimensional depth
model on the same five recovery contrasts and all four known state controls.
There is no choice of spline size, donor subset, effect threshold or alternative
region after seeing outcomes.

The estimand is the conditional confirmation-RNA log-predictive advantage of
two regional populations over one continuous-state population. Each cell's
gene counts follow the same working NB law with archived expected counts and
sizes. Let x be standardized log outside-event discovery depth. The null is
eta=a+b*x+sqrt(v)*Z. The alternative adds delta*K, where Z is standard Gaussian
and K is Bernoulli(pi). The latter shares one prevalence across depths and a
common within-component variance. Both models fit the same trend; delta0 nests
the null. Null/alternative have3/5 fitted parameters rather than separate
locations for each depth stratum.

The identifying restrictions are linear residual-depth response, an adequate
Gaussian state law, transferable reference/count noise, and prevalence
independent of outside-event depth. Nonlinear transcription/depth association,
genotype-associated depth, reference mismatch and RNA-equivalent biological
states remain competing explanations. Independent technical/reference data
could test the depth/noise law; same-cell DNA could test genomic memberships.
This bounded RNA screen supplies neither. Biological n=0; all nine simulations
are exposed regression probes, not biological replication or independent
evaluation. Previously fitted nuisance and adaptive region selection remain
conditioned upon.

Every original contrast must retain a positive numerically resolved advantage
and at least30 discovery responsibility mass in both components. Every state
control must have a negative resolved advantage. Unsupported or unresolved
cases fail the9/9 continuation gate. Failure rejects this unchanged candidate
without a slope/threshold sweep. Passing permits only a separately frozen
full original recovery and fresh whole-workflow false-split/power study. The
5% specimen criterion, exact confidence family, simultaneous power assessment,
all original recovery assertions and frozen artifacts remain unchanged.

## Exact input, inference and evidence chain

[Frozen design](../dev/linear_depth_state_20261111/DESIGN.md),
[new observation model](../dev/linear_depth_state_20261111/model.py),
[fitting/execution driver](../dev/linear_depth_state_20261111/study.py),
[source identities](../dev/linear_depth_state_20261111/source.sha256).

Preparation verifies every existing prepared-array and immutable fold,
confirmation and proposal-lock hash. It retains the original cell/gene order
and selected population. Proposal counts extracted by global row and gene IDs
must exactly equal the archived D matrix. The entire target bins are excluded
from outside-event totals. Log depth is centered and divided by its discovery
RMS; no confirmation observation enters this covariate. A constant unsupported
covariate fails explicitly. The old categorical offsets are undone from both
expected matrices before fitting the new trend, and all offsets are retained.

The field `outside_event_bins` in preparation lists the excluded EVENT bins;
the independent reader labels it `excluded_event_bins`. The implemented mask
uses their complement correctly. This reporting-name limitation is disclosed
without rewriting the already frozen records.

[Input lock](../dev/linear_depth_state_20261111/prepared_v1/INPUTS_LOCKED.json)
records every source, identity, transformation, count axis and covariate scale.
Its separate case arrays contain D and expected confirmation means, but no
confirmation counts. The latter remain in their existing files until each
case's full discovery-model lock is written.

Null fitting uses two fixed variance starts. The alternative uses delta0 and
positive/negative discovery spread crossed with zero/null variance, plus one
native-proposal prevalence start:seven starts. Successful finite L-BFGS-B fits
are local multistart solutions; global optimality is not certified. There is
no variance floor or forced positive copy effect. At variance0, the Gaussian
heat-equation derivative supplies the analytic boundary score. The nested-null
candidate must not lose its feasible likelihood beyond an arithmetic budget.

All three quadrature fits24/48/96 lock before confirmation count decoding.
Conditional prediction is logP(D,C)-logP(D), using one latent state shared by
the two folds. Full gene counts must equal D+C exactly, aligned independently
through stored global IDs and original event genes. This does not establish
thinning independence when the estimated sampling/reference law is inadequate.

Numerical resolution remains twice the largest consecutive-order score
difference plus the operation/magnitude rounding budgets at all three orders.
This is a convergence diagnostic, not a quadrature enclosure, confidence
interval or biological false-split calibration. No convenient epsilon defines
a scientific sign or replaces30-cell support.

## Independent computational prerequisite

The numerical/input prerequisite23206618 completes0:0 on nslab CPU in16seconds.
All nine discovery inputs lock successfully. The
[numerical report](../dev/linear_depth_state_20261111/prepared_v1/NUMERICS.json)
and [independent reference code](../dev/linear_depth_state_20261111/checks.py)
retain the following results:

- Exact Fraction zero-variance NB probabilities and gradients agree to
  1.1102230e-16, with operation/magnitude budget1.1510577e-9.
- Independent adaptive integration of full SciPy NB PMFs and Gaussian
  heat-equation gradients agrees to5.4123372e-16. Adaptive error, explicit
  12-SD tail bounds and rounding give budget1.1562650e-9.
- All five finite-difference derivatives pass step-refinement/rounding budgets.
  Largest observed error4.3409554e-11; its budget6.4628279e-10.
- Nested-null, component exchange, cell order and covariate reparameterization
  invariants pass.

These checks assess computation under the stated model. They cannot establish
that its latent-state law or biological interpretation is correct. No public
package source changed, so the prior521-Python/four-Rust validation is historical;
the full package/GPU suites were not rerun for this separate research model.

## Execution checkpoint

One prerequisite, nine independent case jobs and a dependent summary are
submitted on nslab account/QOS,CPU,one CPU each. Prerequisite8GB/1h;
cases8GB/3h;summary4GB/1h. A separate dependent reader requests8GB/1h.
Every submission batch was announced first. There are no arrays, interactive
QOS, dependency installations, shared-environment changes or running-job edits.

[Prerequisite command](../dev/linear_depth_state_20261111/prepare_submission.json),
[all nine commands](../dev/linear_depth_state_20261111/case_submissions.jsonl),
[dependent summary](../dev/linear_depth_state_20261111/summary_submission.json),
[reader command](../dev/linear_depth_state_20261111/reader_submission.json),
[reader protocol](../dev/linear_depth_state_20261111/READER_DESIGN.md),
[independent reader](../dev/linear_depth_state_20261111/reader.py).

At preparation of this checkpoint,seven of nine cases have completed, two are
still running and summary/reader await dependencies. No complete scientific
gate, recovery gain, false-split improvement or final validation is claimed
from those partial results. All new output directories are exclusive and all
old results remain. Final outcomes and complete accounting will follow below.

## Complete nine-case result:8/9, reject the candidate

All12 jobs complete0:0. The scientific screen fails. No case, support threshold,
numerical-resolution criterion, count law, comparator or outcome is replaced.
The [complete summary](../dev/linear_depth_state_20261111/summary_v1/REPORT.json)
contains every order-specific result and input link.

| Fixed case | Predictive advantage | Numerical diagnostic | Smaller discovery mass | Pass |
|---|---:|---:|---:|---|
| original20261031 | +5.595246 | .00064140 | 217.6517 | Yes |
| original20261107 | +4.821689 | .00091138 | 196.9199 | Yes |
| original20261108 | +5.171751 | .00083630 | 242.8023 | Yes |
| original20261109 | −3.634283 | .00066595 | 96.6425 | **No** |
| original20261110 | +1.856859 | .00059590 | 139.3650 | Yes |
| state202611096001 | −1.415522 | .00189708 | 52.4751 | Yes |
| state202611096002 | −4.073089 | .00082114 | 63.4237 | Yes |
| state202611096004 | −5.704886 | .00125638 | 190.8665 | Yes |
| state202611096006 | −8.944799 | .00102866 | 173.8024 | Yes |

All signs resolve beyond the unchanged diagnostic, and every final discovery
component has mass above30. The required109 contrast fails because the
continuous-state model predicts confirmation better, not because of a tiny
unsupported component or a quadrature sign ambiguity. Its failing sign is
consistent with the earlier fixed/soft-prevalence screens, although model
families and fitted means differ and these scores are not a matched causal
decomposition.

The [109 discovery lock](../dev/linear_depth_state_20261111/case_03_v1/DISCOVERY_MODELS_LOCKED.json)
records null(a,b,v)=(.01316534,−.00087674,.01063341). The fitted alternative has
a=−.03567760,b=−.00041855,delta=.25382451,logit(pi)=−1.43375038,v=0, with
component masses405.3575/96.6425. Its nonvariance gradient infinity norm is
1.01894e-8. These are stored local-optimization diagnostics, not a certified
global solution or a DNA copy contrast. Successful discovery fit and a
nonzero contrast do not establish predictive or genomic support.

This rejects the specific claim that a common linear depth trend is sufficient
to repair the membership/state tradeoff. It does not prove that depth has no
role or that no other justified model can improve recovery. Changing the slope
basis, prevalence constraint or sign threshold after this result is not an
admissible continuation of this candidate. No full-workflow inference or fresh
null/power panel is run for it. Additional computation on this unchanged
correction would not repair its failed continuation gate.

## Independent completion checks

[Reader23206642](../dev/linear_depth_state_20261111/reader_23206642/results/REPORT.json)
independently sums all27 per-cell arrays as exact Fractions, reconstructs every
convergence diagnostic and all9 decisions, and verifies immutable input hashes.
Its full-PMF scalar Gaussian integration probes36 cells selected by first/last
global ID and minimum/maximum outside-event depth, with no favourable-score
selection. The maximum independent conditional-score discrepancy is
3.6970427e-14. Every probe passes its adaptive/tail and operation/magnitude
budget. The normal tail is truncated using the gene-wise maximum PMF and
log-tail inversion, rather than an arbitrary fixed integration range.

These36 checks are bounded numerical references, not a certified enclosure for
all600-cell quadrature calculations or a statistical confidence guarantee.
Original per-case count conservation succeeds exactly after fitting. All27
discovery fits precede their confirmation reads; no all-cell descriptive refit
overwrites the recorded score. The preparation's numerical and input checks
pass separately from the failed scientific gate.

| Job / role | Progress | State | Elapsed | Peak process RSS,KiB |
|---|---|---|---:|---:|
| 23206618 numerics/input | 1/1 ████████ | Complete0:0 | 16s | 107996 |
| 23206621 original1031 | 3/3 ████████ | Complete0:0 | 51s | 467348 |
| 23206622 original107 | 3/3 ████████ | Complete0:0 | 1m09s | 442168 |
| 23206623 original108 | 3/3 ████████ | Complete0:0 | 50s | 448020 |
| 23206624 original109 | 3/3 ████████ | Complete0:0 | 3m23s | 434988 |
| 23206625 original110 | 3/3 ████████ | Complete0:0 | 6m04s | 444408 |
| 23206626 state6001 | 3/3 ████████ | Complete0:0 | 7m34s | 499792 |
| 23206627 state6002 | 3/3 ████████ | Complete0:0 | 8m34s | 509372 |
| 23206628 state6004 | 3/3 ████████ | Complete0:0 | 1m10s | 505288 |
| 23206629 state6006 | 3/3 ████████ | Complete0:0 | 1m47s | 510720 |
| 23206630 complete summary | 9/9 ████████ | Complete0:0 | <1s | 41824 |
| 23206642 independent reader | 36/36 ████████ | Complete0:0 | 6s | 89232 |

Running0,pending0,complete12,executable failures0 in this chain. Scientific
failures1/9. [Scheduler accounting](../dev/linear_depth_state_20261111/accounting_complete.tsv)
and [GNU-time measurements](../dev/linear_depth_state_20261111/resources_complete.json)
retain exact nodes, timings, allocations and memory. Subsecond accounting is
displayed as<1s, not interpreted as no computation. The maximum process peak
is510720KiB; short-job Slurm RSS0 is not used as a memory estimate.

[Integrity checks](../dev/linear_depth_state_20261111/source_integrity_complete.json)
verify all91 previously validated package/native/data/test identities, all
new executed source manifests, and the previous turn's exact arm scorer,
driver fixes and tests. No existing package, benchmark, native or test code
was edited for this experiment. The ledger adds this failed candidate while
retaining all281 previous values. The documentation preserves the earlier
partial checkpoint as history and records complete results here.

The next justified performance intervention must alter information or fitting
where the required109 contrast is lost; this experiment does not support more
categorical-depth or prevalence tuning. The independent numerical agreement
rejects rounding as an explanation for this negative score. No merge ran here,
so this cannot establish or refute a merge correction. Original production
recovery remains2/5, continuous-state false splits14/200 and correlated selected-
interval coverage54/1000. All frozen outputs and scientific criteria remain.
No competitor, biological, spatial, membership-confidence or SOTA improvement
is claimed. The full scientific objective remains active and unachieved.
