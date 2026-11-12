# Confirmation, merging and the estimand of a pooled difference

The scientific goal is to preserve recoverable membership information without
turning uncertain assignments or expression-state differences into supported DNA
clones. Every original recovery assertion and the 5% specimen false-split target
remain fixed. Frozen evidence is retained; these are development simulations,
not independent biological validation. Performance findings are
[Observed pattern, not a guarantee].

## What was changed and what was tested

`SpecimenModel.matched_differences` now accepts explicit nonnegative cell weights
and an optional minimum effective-cell support. Matching balances **weight mass**
within the existing depth strata. The conditional sampling term uses squared
weights; pseudocounts retain the existing weight-unit scaling. Kish support is
computed after matching. Constant weights have exactly their positive cell count,
preventing floating-point rounding from creating a spurious support failure.
The historical unweighted default is retained.

The experimental benchmark uses `abs(2*q-1)` for binary event proposals and
`abs(R_a-R_b)` for merge pairs. These are classification margins, with no fitted
confidence cutoff. It keeps the original initial HMM, posterior threshold,
minimum segment length, depth matching and gene veto. Effective support must
also meet the original `min_cells=30`. Full-count hard-pool initialization and
export scoring are unchanged. Spatial and allele-enabled fits are explicitly
unsupported by this prototype.

The estimand is a confidence-weighted **selected-cell RNA contrast**, not the
latent component contrast. Squared-weight propagation is computationally valid
conditional on fixed weights and the Poisson sampling term. It does not account
for fitted probabilities, NB excess variation, nuisance estimation or selection.
Proposal-half assignment after full-data nuisance fitting is not independent
confirmation. The strongest competing explanation remains expression-state
extremes; the repeated full workflow is needed to test this failure mechanism.

The intervention is isolated in `benchmarks/hardening/membership_weighted_stages.py`.
It records event probabilities before they are discarded by hard assignment;
it uses the existing proposal-half probabilities at merge time. It does not use
truth labels. Profiling hooks are experimental instrumentation, not a new public
clone API. The package change is the tested explicit-weight pooling primitive;
no inference default is promoted by these experiments.

## The five-case causal ablation

All four paths use the same frozen count matrices, positions, truth-free fitting
metadata, parameters, genes, thresholds and scorer. Truth is joined only after
inference. Original successful cases are retained as controls.

| Seed | Hard control | Confirmation only | Merge only | Both |
|---|---|---|---|---|
| 20261031 | pass | pass | pass | pass |
| 20261107 | pass | pass | pass | pass |
| 20261108 | fail | fail | pass | pass |
| 20261109 | fail | fail | fail | fail |
| 20261110 | fail | fail | fail | fail |

For seed 20261108, merge weighting preserves two tumour groups with assignment
accuracy 0.7668, an exported chr12 distinction and normal false calls 8/90.
The hard control has one tumour group and 10/90 normal false calls. Every original
criterion passes in the weighted-merge path. Confirmation weighting alone has no
effect on recovery. This is evidence that hard reassignment before merging can
remove recoverable information; it is not evidence of a reliable general clone
caller. The simpler merge-only intervention therefore remains an explicit
comparator in the error/power experiment.

Seeds 20261109 and 20261110 still fail the initial weighted confirmation HMM for
the chr12 event. Neither reaches the depth-matched confirmation check for that
proposal. Mean classification margins are 0.4158 and 0.5225, respectively. The
failure is not repaired by preserving the same accepted groups at merge time,
because those groups are never accepted.

## Why confidence weights do not fully undo mixing

For a proposed low/high component with rates `theta0, theta1`, a pool with
exposure-weighted high-component probability `p` has expected rate
`(1-p)*theta0 + p*theta1`. A contrast between mixed pools therefore differs from
the latent log-rate contrast. Independent direct averaging and log-sum-exp
calculations on saved proposal fits agree within 5.3e-16 log-rate units.

| Seed | Latent fitted contrast | Expected hard-pool contrast | Expected margin-weighted contrast |
|---|---:|---:|---:|
| 20261108 | 0.2646 | 0.1380 | 0.1789 |
| 20261109 | 0.2443 | 0.1028 | 0.1414 |
| 20261110 | 0.2539 | 0.1393 | 0.1751 |

These are **conditional predictions using fitted probabilities**, not measured
pool purity or independently calibrated membership uncertainty. Weighting reduces
mixing while lowering effective sample size. It does not identify latent profiles
or justify changing the HMM threshold. The current HMM compares the attenuated
pool contrast with whole-step alternatives. A future latent-mixture observation
model would need independent identifiability and uncertainty checks, followed by
the same full-workflow null/power tests. Substituting latent point estimates while
retaining pure-pool standard errors would not be a justified correction.

## Compute validation and provenance

Work is recorded under `dev/confirmation_20261107`. The dirty baseline contains
349 source files captured before editing. Builds and all executable checks run
on compute nodes in a fresh isolated environment. Existing environments and
frozen outputs are not modified. CPU analyses use Intel v4 nodes because earlier
work demonstrated unresolved cross-generation numerical differences.

Job 22855195 passes ten focused baseline Python tests and four Rust tests.
Candidate job 22855685 passes **319 Python tests, one GPU-only skip, and four
Rust tests**, including the original primary recovery regression. Its peak RSS
is approximately 2.13 GiB. New checks independently verify fractional mass
matching, squared-weight variance, exact legacy/unit-weight equality,
weight-scale invariance, effective support and zero-information abstention.
Numerical comparisons use machine-precision tolerances tied to these small
algebraic calculations, not a scientific decision threshold.

Recovery jobs 22855713–22855716 and independent expectation job 22856159 complete
successfully. `candidate_v1`, `candidate_v2`, `candidate_v3` and subsequent audit
snapshots preserve exactly executed source. Later snapshots add trace identity
or summary/audit logic; the tested installed package is unchanged. Original
scientific tests and simulation generators are preserved.

## The proposal prediction also needs checking

A bounded follow-up uses the saved event-proposal exposure, NB size, component
rates and probabilities to predict the corresponding confirmation-region counts.
It integrates each cell's shared latent component when computing total variance:
the within-component NB variance plus the variance between component means.
Independent SciPy NB moments reproduce this calculation. This is a conditional
model check, not a calibrated interval or a multiplicity-adjusted hypothesis test.

| Seed | Predicted margin-pool contrast | Observed confirmation contrast |
|---|---:|---:|
| 20261108 | 0.1789 | 0.1445 |
| 20261109 | 0.1414 | 0.0984 |
| 20261110 | 0.1751 | 0.1001 |

These regional totals precede the HMM's programme correction, centering and
segmentation, so they are not the final HMM contrast. Job 22857115 records the
actual counts, exposures, conditional variances and post-fit truth composition.
For example, in seed 20261110 the fitted high-component fraction in the weighted
high group is 0.849, while the simulated-truth fraction is 0.748. The proposal
probabilities cannot simply be treated as known genomic membership fractions.
Selection, nuisance estimation and ordinary sampling all remain possible causes;
these three selected cases do not estimate their population frequency.

The smallest next correction must therefore model the observable mixture and
propagate probability/nuisance uncertainty, or abstain. Merely dividing the
observed contrast by fitted attenuation, hardening the assignments further, or
lowering the HMM threshold is not supported. An independently generated
confirmation-prediction check with known and estimated nuisance quantities would
distinguish an implementation problem from estimation and selection error.
A model-based success would still require biological validation at the claimed
clone resolution.

## Complete binned-workflow false-split and power result

All three variants (hard, merge-only, both) complete the same 998 fixed simulations:
199 draws from each of two calibration laws and 200 evaluation draws in each of
three conditions. Across variants this is 2,994 fits. Counts, known reference,
dispersion, latent state and simulated truth are paired exactly. The workflow
includes initial groups, proposal search, confirmation, merging, splitting and
stopping. The score is the original maximum supported final-group difference,
not the earlier residual-event observer score. Public raw-gene preprocessing,
reference estimation, gene veto and alleles are outside this binned experiment.

The rank level remains 0.05. Exact one-sided binomial bounds use Bonferroni
0.05/27 across all variants, decision rules and conditions. Adding the simpler
merge-only comparator enlarged the family before null results were summarized;
this tightened the criterion. Original protocol and addendum are retained.
Bounds condition on the realized shared calibration banks. The replication unit
is an independently simulated specimen, not a cell or patient.

All three variants have the following marginal counts and bounds:

| Decision rule | Matched-NB false distinctions / 200 | State-null false distinctions / 200 | Alternative detections / 200 |
|---|---:|---:|---:|
| Raw multigroup output | 0; upper 3.10% | 190; upper 98.37% | 200; lower 96.90% |
| NB-calibrated final-group gate | 0; upper 3.10% | 190; upper 98.37% | 200; lower 96.90% |
| State-calibrated final-group gate | 0; upper 3.10% | 4; upper **6.81%** | 65; lower **23.28%** |

**The 5% error criterion fails.** Equal marginal counts do not mean every decision
is identical: each weighted variant gains one alternative detection and loses
one compared with hard pooling. All 18 paired exact discordance tests have
Holm-adjusted p=1. There is no demonstrated improvement in this error/power
tradeoff, and no equivalence claim is made.

The raw alternative is also heavily over-split: exactly two eligible tumour
groups occur in only 1/200 hard-control specimens and 0/200 weighted specimens.
Median raw membership ARI is 0.07646 (hard), 0.07554 (merge), and 0.07590 (both).
Global detection is therefore not correct clone recovery. The calibrated gate
does not validate all retained pairs or any individual membership. The matched
NB success does not establish robustness to correlated expression state.

The merge change recovers one consequential development case without losing an
original success, but neither weighted variant is promoted to an inference
default. The explicit-weight pooling primitive and reproducible experimental
paths remain available; two confirmation failures and the error-control failure
remain open. No biological specimen or suitable unopened validation panel was
added, so the proposed biological guarantee remains unsupported.

## Independent checks and delivery

Summary job 22856898 checks all paired inputs and complete shard populations.
Audit job 22856899 independently recounts decisions and inverts binomial tails;
its bounds differ from the beta-quantile calculation by at most 1.81e-16.
It verifies all five hard-control posterior arrays and scientific tables exactly
against the prior frozen panel. It also checks 349 initial source files, 349
previous-final source files, the original criteria/generators, 25 panel-manifest
input records, nine selected historical result records, eight prior build
artifacts and 22 prior installed package files. This is a bounded preservation
audit, not an exhaustive inventory of historical raw data.

Identity job 22858094 independently verifies **all 998 prior hard-control input
records, probability arrays, assignments and original final-group scores exactly**.
The previous conditional event-observer score is deliberately not used as the
new final-group score. All 27 jobs in this turn complete successfully; scientific
failures are retained in their reports rather than hidden as execution successes.

Changed source is listed in `dev/confirmation_20261107/change_manifest.json` and
separated from the pre-existing dirty tree in `changes.patch`. `final_source`
preserves the delivered source; `jobs_final.tsv` records accounts, resources,
runtimes, exit codes and dependencies. The tested noneditable wheel SHA256 is
`93205a319db1c1000148bf82efa8dd4f14a7c60f492d73e1a48c485123f833e9`.
