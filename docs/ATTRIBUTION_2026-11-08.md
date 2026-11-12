# Membership and profile-model attribution

The goal remains accurate RNA-derived dosage profiles and reliable CNA-defined
clones, with superiority established by appropriate frozen comparisons rather
than asserted from design features. This round keeps profile/bulk and clone
reliability in view. No inference default, threshold, original scientific test,
former evaluation result or prior execution snapshot changes.

## Clone confirmation: two sources of error

The 18-fit factorial diagnostic crosses fitted soft, hard and known simulated
membership with fitted and oracle exposures for the three locked failures.
The simulator must reproduce every frozen CSR count and cell ID exactly. True
gene means are aggregated using the executed gene/bin mapping. Oracle exposure
removes the known regional effect from each cell's true mean; it includes
population information and is explicitly nondeployable. Truth never enters
HaploDose inference. NB size remains fixed, isolating the mean correction without
pretending that unequal-gene aggregation or estimated thinning becomes exact.

| Seed | Fitted exposure / hard | Fitted / soft | Fitted / truth | Oracle / hard | Oracle / soft | Oracle / truth |
|---|---:|---:|---:|---:|---:|---:|
| 20261108 | .1093 | .2411 | .2417 | .1127 | .2489 | .2531 |
| 20261109 | .0742 | .1720 | .2194 | .0794 | .1752 | .2250 |
| 20261110 | .0638 | .1872 | .2190 | .0685 | .2024 | .2316 |

The known simulated log effect is .2433. These are conditional likelihood fits,
not end-to-end recovery or calibrated interval estimates. In the two unresolved
cases, replacing fitted memberships by truth changes the contrast more than
replacing fitted exposures alone. Hard labels retain substantial attenuation
even with oracle exposure. This supports membership loss as a consequential
mechanism while retaining nuisance estimation, thinning and aggregation as
additional sources of error [Observed pattern, not a guarantee].

The likelihood is checked with independent BFGS and Nelder-Mead optimizations;
hard-membership estimates also have independent scalar score roots. Precision
loss flags are retained. Agreement between optimizers cannot establish global
optimality, statistical coverage or biological correctness. Conditional SEs omit
q/nuisance estimation and selection and are not promoted to calibrated intervals.

## Profile/bulk: independently selected segmentation boundaries

**Correction:** the TCGA shared-penalty result below used incompatible log units
and is superseded by [the unit-corrected replay](MIXTURE_CONFIRMATION_2026-11-08.md).
The original numbers remain here as historical evidence. The other four cases
are unaffected; the TCGA claim that naive segmentation wins is withdrawn.

A shared penalized least-squares operator is applied independently to naive and
HaploDose RNA profiles in all five previously locked cases. Gaussian RNA log
ratios are a working loss, not a raw-count uncertainty model. The common noise
scale comes from naive RNA first differences, without DNA or HaploDose input.
The fixed description-length penalty is `3*sigma^2*log(n)` per segment: mean
parameter cost plus breakpoint-location cost, up to constants. Serial RNA
correlation can violate that working model; no error-control claim follows.
The existing minimum segment length is retained. Chromosome arms, missing RNA
support and nonconsecutive bins break the fit. Short blocks receive a descriptive
constant. DNA enters scoring only after boundaries are selected.

An exhaustive enumeration of every partition of a seven-value example verifies
the dynamic-programming objective at three minimum lengths. This control
addresses the earlier asymmetric use of HaploDose-selected boundaries for both
methods. It is still a small development diagnostic, not a population comparison.

| Locked case | Naive, own boundaries | HD ratio, own boundaries |
|---|---:|---:|
| 2313287_STOMACH | .7526 | .6650 |
| 42MGBA_CENTRAL_NERVOUS_SYSTEM | .8130 | .7961 |
| U-CAN U0001 | .7669 | .7906 |
| U-CAN U0008 | .5910 | .6355 |
| TCGA-05-4389-01 | .8038 | .7847 |

Values are Spearman agreement on identical callable DNA support. HaploDose's
ratio advantage survives this operator control in the two U-CAN tumours, but
not in the two Kinker lines or transformed-input TCGA case [Observed pattern,
not a guarantee]. The TCGA case cannot validate raw-count uncertainty. There
are no biological-unit p-values or general superiority claims from these five
selected examples. The simple baseline remains mandatory.

## What remains required

Computational correctness: preserve the verified count/coordinate/numerical
repairs; continue to test any new observation model against independent numeric
references. This round introduces diagnostic code only.

Clone reliability: test an observable-mixture confirmation model that retains
uncertain memberships and accounts for nuisance fitting. It must meet the same
original recovery assertions and complete-workflow null/power criteria. The last
state-calibrated false-split upper bound remains 6.81%, exceeding 5%; neither this
oracle experiment nor a small conditional SE changes that result.

Profile/bulk accuracy: isolate reference, programme and depth corrections against
these independently segmented simple baselines. Keep mode-specific improvements
only when paired evidence supports them. Better smoothing alone is not evidence
of a better genomic signal model.

Uncertainty and biological validation: routine selected-interval coverage remains
unresolved; an oracle-covariance procedure does not supply missing reference or
covariance knowledge. Same-cell or suitable resolution-matched DNA and adequately
replicated unopened specimens remain necessary for biological clone claims.

Spatial assignments and competitor superiority: spatial mixing/boundary controls
and independent validation remain open. Existing inspected cohorts are regression
only. A SOTA claim requires a declared competitor family, task-specific endpoints,
coverage/failure rules and a frozen candidate evaluated independently. No result
in this round establishes that claim or supports universal superiority.

## Execution records

All numerical work runs on Intel v4 compute nodes under innovation, using the
previously verified wheel read-only. `dev/attribution_20261108/baseline` preserves
358 source files and the dirty Git state. New output directories and immutable
source snapshots bind each executed script and scorer. Job 22961013 completes
the clone attribution; 22961071 completes the profile operator control. The
original 319-pass/one-skip Python and four-pass Rust result remains the latest
package validation; a new full suite was not run because package code is unchanged.

Audit job 22961139 verifies both 358-file frozen source inventories and all five
frozen profile inputs. The nine BFGS precision-loss warnings remain recorded;
all independent Nelder-Mead runs succeed. Maximum objective disagreement is
1.82e-12 and maximum hard-membership parameter disagreement with independent
scalar score roots is 4.21e-12. These discrepancies are negligible relative to
the observed contrast differences but do not remove model/selection uncertainty.
All three jobs complete successfully. Peak resident memory and runtime are
recorded in each analysis job's `resources.txt`; scheduler RSS is unavailable
for these short jobs and must not be interpreted as zero memory use.
