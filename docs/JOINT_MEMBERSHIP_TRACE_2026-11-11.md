# Where the low-depth joint model loses profile information

The full genomic and competitor objectives remain active. The preceding joint
state/profile model improved the simulated state-free target at high depth but
failed seven low-depth comparisons. This continuation isolates whether hard
component selection, profile estimation or uncertain memberships causes those
failures before choosing another estimator. It preserves the full performance,
false-split, power and uncertainty criteria.

[Previous full results](JOINT_STATE_PROFILE_2026-11-11.md),
[fixed diagnostic design](../dev/joint_membership_trace_20261111/DESIGN.md),
[immutable source identities](../dev/joint_membership_trace_20261111/SNAPSHOT.json).

## Measurement and identification contract

Use all 288 original exposed multinomial RNA cases: three scenarios, two depths
and 48 generating seed units per scenario, paired across depths. Biological n is
zero. Read the frozen inputs and fitted one/two-component models. No new fit,
feature selection, fold, threshold, state factor or regularisation is introduced.
The existing model is native to these conditional RNA counts; it is not the
heterogeneous gene-NB observation law used by the original false-split panel.

The RNA estimand remains expected regional log odds. The state-free estimand
sets the fitted continuous expression coordinate to zero and averages over
discrete-component responsibilities. It is a model counterfactual, not measured
DNA. Gaussian state independent of population, adequate loading and reference
correspondence remain identifying assumptions. Discrete transcriptional state
can produce the same observations as the proposed genomic profile. Only suitable
independent DNA can distinguish those explanations biologically.

## Decisive decomposition and expected outcomes

The new count-only stage exports the already fitted two-component estimates
even when the original predictive selector chose one. It must reproduce every
one of the 9,216 saved selected-model arrays exactly. All 288 new prediction
bundles lock before truth-based calculations.

The independent reader checks complete multinomial probabilities and posterior
arrays, then compares the original selected estimator with four diagnostic
actions: forced two-component estimates; correct memberships with fitted
profiles; correct profiles with fitted memberships; and both oracles. One
global label permutation uses all-cell MAP agreement. Exact responsibility
ties cannot choose its orientation. Non-tied alignment must be invariant to
reversing component labels. Both oracles together must give exactly zero error.

Improvement from forced two implicates hard selection. Improvement from correct
memberships implicates membership estimation; improvement from correct profiles
implicates profile estimation. Requiring both oracles indicates interacting
errors. These are optimistic diagnostics, not deployable methods. They cannot
establish an information limit if the oracle comparison is inconclusive.

Retain every control, whole/minority endpoint, homogeneous forced-two failure,
high-depth result and numerical boundary overlap. Summaries are paired
descriptive effects over the fixed 48 seed units, with no additional nominal
population tests. Any proposed correction must subsequently pass the original
32 RNA and simultaneous 64 both-target criteria. The public minimum stays 30;
the 24-cell toy minority remains ineligible for supported genomic clone calls.

## Compute and preservation

Trace **23218248** and dependent reader **23218249** each request one CPU,
8 GB and one hour on nslab account/QOS, CPU AMD-v6. Existing isolated Conda and
runtime003 only, one thread, source hashes before/after, neutral directories,
process resources and scheduler accounting. No arrays, new dependencies, public
source changes or running-job modifications. Both completed `0:0`: trace runtime
45.38 seconds, peak process RSS 89,980 KiB; reader 13.72 seconds, 103,604 KiB.

[Exact submissions](../dev/joint_membership_trace_20261111/submissions.json),
[cached-estimate export](../dev/joint_membership_trace_20261111/trace.py),
[independent decomposition](../dev/joint_membership_trace_20261111/analyze.py).

## Completed result: forcing two components is insufficient

All 9,216 selected-model arrays replay exactly. Independent full-PMF checks
verify 1,728 forced-two posterior arrays, with maximum discrepancy 4.19e-12.
All 384 heterogeneous double-oracle checks have exactly zero error. Non-tied
oracle alignments satisfy component-label invariance. No new fit was made.

The table gives mean state-free MSE over the same 48 generating seed units at
depth 300. Correct memberships/profiles are diagnostic oracles with frozen
parameters, not methods deployable without truth.

| Case | Action | Whole MSE | Minority MSE |
|---|---|---:|---:|
| Two profiles | Original selected model | 0.00407038 | 0.08293478 |
| Two profiles | Forced two | 0.00553435 | 0.06144888 |
| Two profiles | Correct memberships, fitted profiles | 0.00714014 | 0.06263968 |
| Two profiles | Correct profiles, fitted memberships | 0.00384822 | 0.06703954 |
| Two profiles plus state | Original selected model | 0.00429165 | 0.08490044 |
| Two profiles plus state | Forced two | 0.00685239 | 0.06417825 |
| Two profiles plus state | Correct memberships, fitted profiles | 0.00695856 | 0.06495383 |
| Two profiles plus state | Correct profiles, fitted memberships | 0.00465330 | 0.04945489 |

Forced two lowers minority error but increases whole-profile error in both
conditions. Correcting memberships alone does not fix the estimated profiles;
correcting profiles alone does not resolve both endpoints. [Inference] The
paired ablations implicate interacting profile-estimation and membership errors,
rather than a simple component-admission mistake. They do not establish an
information limit or prove that any particular penalty will work.

Homogeneous forced-two error also increases: state-free whole MSE is 0.00383361
versus selected 0.00092882 at depth 300, and 0.00033819 versus 0.00007904 at
depth 3,000. High-depth heterogeneous forced-two estimates reproduce the selected
ones, as expected because every selected model already has two components.
All original controls, negative results and numerical diagnostics remain.

[Complete independent result](../dev/joint_membership_trace_20261111/reader_23218249/results/REPORT.json),
[trace resources](../dev/joint_membership_trace_20261111/trace_23218248/resources.txt),
[reader resources](../dev/joint_membership_trace_20261111/reader_23218249/resources.txt).

The next correction tests the smallest ordered pooling constraint on the same
count likelihood, with both pairing origins fixed in advance and matched simple
controls. It is motivated by estimation error, not a lowered admission threshold.
[Exact adjacent-pooling design](../dev/adjacent_joint_profile_20261111/DESIGN.md).

The 333 prior ledger values and 91 protected production/data/test identities
remain unchanged. Research recovery 5/5, production recovery 2/5, frozen state
errors 14/200 and selected correlated coverage 54/1,000 remain separate from
this profile diagnostic. Biological replication, independent genomic validation
and complete comparator evidence remain missing.
