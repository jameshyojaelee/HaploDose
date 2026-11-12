# Joint continuous state and discrete RNA profiles

The scientific goal is accurate relative dosage profiles and supported candidate
populations, with genomic interpretation established by independent measurements
at the claimed resolution. The preceding marginal RNA predictor reduced RNA
error while increasing state-free response error. This study tests the smallest
model extension addressing that failure: one shared continuous expression factor
and a discrete population profile fitted under one count likelihood. It is an
isolated research implementation, not a production release or DNA validation.

[Frozen design](../dev/joint_state_profile_20261111/DESIGN.md),
[source identities](../dev/joint_state_profile_20261111/SNAPSHOT.json),
[preceding target falsification](PREDICTIVE_PRIOR_CHOICE_2026-11-11.md).

## Measurement and identification contract

Conditional on each retained RNA total, regional counts follow a multinomial
working model with probabilities `softmax(mu + z * delta + t * loading)`.
There are one or two discrete components; `t` is standard Gaussian and independent
of `z`. All parameters share the marginal count likelihood. A fixed final-coordinate
gauge removes the softmax constant ambiguity. This law matches the small profile
generator, not the heterogeneous gene-level negative-binomial false-split panel.
The model must not be transferred across those observation laws without adaptation.

The RNA estimand is posterior expected regional log odds. The separate state-free
estimand sets `t=0` and averages over conditional component responsibilities.
State-free response is a model counterfactual, not DNA copy number. Zero-mean
state, independence of state and component, and an adequate shared loading are
identifying restrictions. A discrete transcriptional population can produce the
same observations as a genomic alteration. Corresponding DNA must distinguish
those explanations; genomic membership requires corresponding-cell evidence.

## Fixed decisive experiment

Six exposed inputs use the first lexical identifier of each original scenario
at depths 300 and 3,000: homogeneous, two CNA-shaped profiles, and two profiles
with continuous state. The counts, folds and references come from locked input
archives. Two additional matched negative controls have continuous state and one
CNA-shaped population. Their seed is fixed before fitting. Biological n is zero.
The eight cases investigate mechanisms; they cannot estimate population error or
power and do not replace the frozen false-split experiment.

Fit one/two components with/without the continuous factor. Select the component
count separately for each route by the sum of both directions of held-out-cell
count prediction. Then refit all cells descriptively, retaining the frozen choice
and predictive statistics. Initialisations, optimiser, folds and quadrature
orders are fixed; failed starts remain recorded. Neither predictive selection nor
all-cell responsibilities are independent genomic confirmation.

All eight RNA bundles must be locked and verified before the reader reconstructs
truth. Report both targets, every old control, minority-region errors and the
matched no-factor ablation. Report responsibility masses and MAP counts separately.
The true minority has 24 cells, below the unchanged public support minimum of 30;
these latent profile components are not supported genomic-clone calls.

The initial mechanism gate requires one component for both negative controls at
both depths. The high-depth state-plus-profile case must select two components
and improve whole state-free error over pseudobulk and the matched no-factor
route, while preserving minority accuracy against the old empirical estimator.
Differences must clear the numerical diagnostic. Failure rejects this fixed
candidate without adding factors, cutoffs or favourable starts. A pass permits
a separately frozen full regression; it does not waive any original performance,
false-split, power, uncertainty, comparator or biological criteria.

## Independent checks and compute record

Qualification job **23217686** completed `0:0` on `ne1dc6-003`, nslab account/QOS,
one CPU, 4 GB, one-hour limit. Runtime was 8.86 seconds and peak process RSS was
59,808 KiB. Independent complete multinomial polynomials, complex-step derivatives
with a halved-step check, both posterior estimands, component-label and factor-sign
invariants, extreme log odds, zero-library rejection and feasible nested objectives
passed. Maximum gradient discrepancy was 2.22e-16. Bounds derive from arithmetic
operation/magnitude counts, not a universal epsilon. Production sources are
unchanged; the historical whole-package suite was not rerun.

[Numerical report](../dev/joint_state_profile_20261111/numerical_23217686/results/REPORT.json),
[measured resources](../dev/joint_state_profile_20261111/numerical_23217686/resources.txt),
[numerical submission](../dev/joint_state_profile_20261111/jobs_numerical.json).

Eight separate case jobs **23217727–23217734** request one CPU, 8 GB and six hours
each on nslab CPU AMD-v6. Independent reader **23217735** requests one CPU, 8 GB
and two hours after successful completion of all eight cases. No arrays, new
dependencies, shared installs, running-job changes or frozen-output writes.
Source hashes are checked before and after execution. The reader independently
reconstructs full multinomial probabilities, exact score sums, high-order
posteriors and adaptive integrals at consequential cells. Convergence across
quadrature orders determines numerical resolution; unresolved comparisons cannot
support a scientific decision. All ten jobs completed `0:0`. Case runtimes were
9.04–13.65 seconds, process peaks 95,848–108,768 KiB. The reader took 1.98 seconds,
92,980 KiB. The scheduler reports zero RSS for short batch steps, so these memory
figures come from each process's saved `time -v` record.

[Scientific submissions](../dev/joint_state_profile_20261111/jobs_scientific.json),
[case catalogue](../dev/joint_state_profile_20261111/case_jobs.json),
[research model](../dev/joint_state_profile_20261111/model.py),
[independent reader](../dev/joint_state_profile_20261111/reader.py).

## Completed mechanism result and limitations

The prespecified mechanism gate passes. Both homogeneous cases and both
state-only cases select one component. The high-depth state-plus-profile case
selects two and improves state-free whole error beyond its numerical diagnostic
against pseudobulk and the matched no-factor route. Its minority error also
improves against the frozen empirical control. These are exposed simulation
outcomes, not error-rate or biological estimates.

| Scenario | Depth | Joint components | State-free whole MSE, joint | Pseudobulk | Matched no-factor | Matched assignment, joint / no-factor |
|---|---:|---:|---:|---:|---:|---:|
| Homogeneous | 300 | 1 | 0.00093781 | 0.00094121 | 0.00093400 | 120 / 120 |
| Homogeneous | 3,000 | 1 | 0.00009653 | 0.00009770 | 0.00009685 | 120 / 120 |
| Two profiles | 300 | 1 | 0.00394449 | 0.00398857 | 0.00393333 | 96 / 96 |
| Two profiles | 3,000 | 2 | 0.00021310 | 0.00355823 | 0.00021334 | 120 / 120 |
| Two profiles plus state | 300 | 1 | 0.00386250 | 0.00398463 | 0.00992555 | 96 / 63 |
| Two profiles plus state | 3,000 | 2 | 0.00027033 | 0.00363314 | 0.01085885 | 119 / 66 |
| State only | 300 | 1 | 0.00109375 | 0.00131970 | 0.01023955 | 120 / 84 |
| State only | 3,000 | 1 | 0.00017243 | 0.00032425 | 0.00816267 | 120 / 64 |

All assignment counts use the same 120 cells and one global label permutation.
They are descriptive all-cell RNA fits. The high-depth mixed case's minority
mass is 23.47, with 23 MAP cells; no supported genomic clone is emitted. At low
depth both true-profile cases collapse to one component, losing the minority.
Homogeneous RNA MSE worsens versus pseudobulk: 0.003565 versus 0.000941 at depth
300 and 0.000264 versus 0.0000977 at depth 3,000. These negatives remain visible.

In the high-depth mixed case, joint state-free whole MSE is 0.00027033 with a
conservative numerical MSE allowance of 0.00221453, versus pseudobulk 0.00363314
and no-factor 0.01085885. Minority state-free MSE is 0.00222899 with allowance
0.01157714 versus empirical 0.04277210. Thus the prespecified differences clear
the diagnostic despite appreciable order sensitivity. This allowance is a
quadrature/optimisation diagnostic, not statistical uncertainty. The largest
independent posterior discrepancy is 3.74e-12; the largest adaptive integral
log-probability discrepancy is 0.000108 and remains within the saved convergence
diagnostic. The full reader retains all cells, controls, sums and checks.

[Complete independently checked outcomes](../dev/joint_state_profile_20261111/reader_23217735/results/REPORT.json).

The mechanism pass permits the separately frozen 288-case regression, using the
byte-identical model and unchanged fitting policy. Eight jobs 23217873–23217880
each request one nslab CPU, 8 GB and three hours for 36 fixed cases; independent
reader 23217881 requests one CPU, 8 GB and two hours. Every original 32 RNA test
remains; add the corresponding state-free tests with simultaneous assessment
across all 64. No new factor, threshold, optimisation start or favourable seed.
All nine full-regression jobs completed `0:0`.
[Full design](../dev/joint_state_profile_full_20261111/DESIGN.md),
[full source snapshot](../dev/joint_state_profile_full_20261111/SNAPSHOT.json),
[exact submissions](../dev/joint_state_profile_full_20261111/submissions.json).

## Full 288-case result: substantial high-depth gain, incomplete overall gate

The unchanged candidate passes **30/32 original RNA comparisons** and **57/64
simultaneous RNA/state-free comparisons**. The full gate fails. All 32 high-depth
comparisons pass under the jointly protected 64-test family. This is a qualified
simulation result at the prespecified depth, not permission to discard low-depth
cases or claim general superiority. Every original criterion stays in force.

Each table entry is the mean over the same 48 generating seed units, paired
across depths within scenario. Cells are not independent biological replicates.

| Scenario | Depth | RNA MSE, joint | RNA MSE, pseudobulk | State-free MSE, joint | State-free MSE, pseudobulk | Joint selects two / 48 |
|---|---:|---:|---:|---:|---:|---:|
| Homogeneous | 300 | 0.00282235 | 0.00091601 | 0.00092882 | 0.00091601 | 0 |
| Homogeneous | 3,000 | 0.00026764 | 0.00007945 | 0.00007904 | 0.00007945 | 0 |
| Two profiles | 300 | 0.00431896 | 0.00413259 | 0.00407038 | 0.00413259 | 0 |
| Two profiles | 3,000 | 0.00034060 | 0.00348458 | 0.00016762 | 0.00348458 | 48 |
| Two profiles plus state | 300 | 0.00725041 | 0.01619652 | 0.00429165 | 0.00440117 | 3 |
| Two profiles plus state | 3,000 | 0.00056371 | 0.01554463 | 0.00036672 | 0.00375644 | 48 |

The high-depth state-plus-profile result corrects the preceding target failure:
the old predictive prior's state-free error was 0.01066673, versus this joint
model's 0.00036672 on identical inputs. Matched no-factor error is 0.01056626.
The joint model therefore improves the simulated state-free target as well as
RNA prediction. Matched all-cell assignments increase from 3,475/5,760 to
5,758/5,760, aggregated descriptively over 48 simulations. Corresponding DNA
and biological validation remain absent. All high-depth minorities remain
below the unchanged 30-cell minimum; no supported genomic clone calls result.

All seven failed simultaneous comparisons occur at depth 300:

| Target | Scenario | Endpoint | Baseline | Numerically resolved wins / 48 |
|---|---|---|---|---:|
| RNA | Two profiles | Whole | Pseudobulk | 24 |
| RNA | Two profiles plus state | Minority | Beta expected log odds | 20 |
| State-free | Two profiles | Minority | Global shrinkage | 16 |
| State-free | Two profiles | Minority | Beta expected log odds | 0 |
| State-free | Two profiles plus state | Minority | Pseudobulk | 3 |
| State-free | Two profiles plus state | Minority | Global shrinkage | 2 |
| State-free | Two profiles plus state | Minority | Beta expected log odds | 7 |

The one-component state-free action removes cell-specific variation, so it
cannot preserve an unresolved minority. At depth 300 it selects one component in
all 48 clean two-profile cases and 45/48 state-plus-profile cases. Homogeneous
RNA error also increases. That generator has no continuous factor; inferred
cell-specific factor variation there comes from sampling. These findings reject
this fixed method as a general replacement.
They do not justify changing support, win-rate, error or uncertainty criteria.
No factor, prior, threshold, fold or start sweep follows this failure.

The independent full reader verifies all 288 input/truth identities, 2,304
archived control rows, 5,184 exact score sums, 1,152 paired held-out probability
checks, 1,728 posterior-array checks and 1,710 adaptive probability probes.
Maximum independent posterior discrepancy is 4.21e-12. Maximum adaptive
log-probability discrepancy is 0.002777, within the saved order-convergence
diagnostic; this diagnostic is not a statistical confidence interval. All 64
Beta-quantile decisions agree with exact rational binomial-tail decisions.

Fit blocks take 4m03s–7m08s with process peak RSS 76,080–89,768 KiB. Reader
23217881 takes 1m17.15s and 104,132 KiB. Combined with the initial qualification
and mechanism panel, all **19 jobs** complete `0:0` on nslab CPU. No new public
package tests are claimed; public sources remain unchanged.

[Complete independently checked full result](../dev/joint_state_profile_full_20261111/reader_23217881/results/REPORT.json),
[reader resources](../dev/joint_state_profile_full_20261111/reader_23217881/resources.txt),
[existing Conda environment inventory](../dev/joint_state_profile_full_20261111/environment_inventory.json).

All 329 pre-existing claim-ledger values and 91 protected production/data/test
identities are unchanged. Both new execution snapshots retain all 22 files
exactly. The earlier inefficient reader 23213185 finished naturally after
1h29m47s; it was not interrupted, resubmitted or modified. Its duplicate work
does not change these scientific results.
[Slurm accounting](../dev/joint_state_profile_full_20261111/accounting_complete.tsv).

## Relationship to the scientific review

This addresses findings 1–5, 7, 9 and 12 by testing the actual target failure and
one common likelihood. It does not resolve biological identification, measured
DNA accuracy, complete-workflow false-split control or selected-interval coverage.
The review's direct code corrections and their executed 521-Python/4-Rust-test
validation are recorded separately; they are not new checks in this study.
[Review response](SCIENTIFIC_REVIEW_RESPONSE_2026-11-10.md).

Research recovery remains 5/5, production recovery 2/5, frozen continuous-state
errors 14/200, and correlated selected-interval coverage 54/1,000. Missing
biological controls and CopyKAT/SCEVAN family outputs remain evidence gaps.
Neither more compute nor a successful mechanism test establishes those claims.
