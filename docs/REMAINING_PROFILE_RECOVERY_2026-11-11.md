# Remaining profile recovery: selection versus collapsed fitted profiles

The scientific goal is accurate genomic profiles and supported clone membership,
with biological attribution established by independent measurements. This
continuation returns to the unresolved low-depth profile/membership failures
without tuning the consumed 50-model DNA-response panel. The preceding goal
turn made progress: it completed that fixed test, rejected promotion and
established current exposure. The full genomic/clone/spatial/SOTA goal remains
active and unachieved.

## Exact LOO candidate and decisive qualification

The corrected ordered joint model passes all high-depth comparisons but fails
seven low-depth comparisons in its full 32/64 family. Its two-fold selector
fits population parameters on 60 cells before a 120-cell descriptive refit.
The hypothesis was that this training-size loss suppresses minority detection.
The bounded candidate implements exact leave-one-cell-out fitting on 119 cells,
with the existing count model, path prior, projection, starts, optimizer and
tolerances. It does not add a fold-count sweep or change a scientific cutoff.

The estimands are conditional expected RNA log odds and the working state-zero
RNA response. Counts follow the existing multinomial/Gaussian-factor mixture.
Independent state/component structure remains an identifying assumption.
Continuous transcription or regional regulation remains a competing explanation
for a CNA-shaped profile. Same-cell DNA or another defensible independent
membership measurement would be needed for genomic attribution. These exposed
simulation cases have biological n=0.

Qualification **23226980** completes 0:0 on nslab CPU in 26 scheduler seconds,
with peak process RSS **97,028 KiB**, requesting 1 CPU / 4 GB / 1 h. It verifies
96 individual/subset/batch score comparisons with independent complete SciPy
multinomial PMFs, exact old batch score arrays, a single cell with zero-count
features, and rejection of six invalid inputs. It reruns all 12 excluded-cell
training fits after perturbing only the excluded cell: every fitted record
remains exactly identical. The training validator and production source remain
unchanged. This is computation and leakage qualification, not calibrated
genomic confirmation.

The new query score helper accepts a legitimate individual cell with unobserved
features; the old research training validator requires two cells and positive
pooled support. That old training requirement is not imposed on held-out queries.
The helper computes the same marginal count score, with unchanged factor
integration and operation-derived arithmetic budgets.

[Frozen LOO contract](../dev/loo_joint_profile_20261111/DESIGN.md),
[exact exclusion/score implementation](../dev/loo_joint_profile_20261111/selection.py),
[qualification and all cached traces](../dev/loo_joint_profile_20261111/qualification_23226980/results/REPORT.json),
[exact submission](../dev/loo_joint_profile_20261111/qualification_submission.json).

## Why the expensive LOO mechanism panel was not submitted

The same qualification checks all 290 cached cases. In each first lexical
low-depth heterogeneous case, the **full 120-cell joint K2 fit already has
identical components**, with responsibility masses approximately 60/60 and
no resolved training count-likelihood improvement. Its gains are 9.44e-12
without added state and 5.57e-12 with state, consistent with arithmetic noise.
The two-fold predictive gains are likewise unresolved. No genuine membership
information remains in those cached full fits for a new selector to recover.

This falsifies the sufficiency of a selector-only correction **for these
ordered-model fits**. It does not prove an information limit for the measurements,
rule out another observation model or certify global optimization. The known-
parameter count oracle still passes, and supplying known parameters changes
the available information. Do not reinterpret that oracle as a deployed method.

The LOO writer and launcher were prepared but **not executed**. No eight-case
or 290-case LOO result exists. No mechanism pass or performance gain is claimed.
Its additional low-depth criterion would have required improved state-free
minority error as well as the original mechanism gate; it was not weakened.

[Prepared, unexecuted writer](../dev/loo_joint_profile_20261111/study.py),
[prepared, unsubmitted launcher](../dev/loo_joint_profile_20261111/case.sbatch).

## Next decisive distinction: remaining search miss or unsupported complexity

A better feasible point under the existing objective would prove a remaining
optimization failure; absence of such a point would not prove global optimality
or identify the prior as the cause. The next job compares all 5,220 cached
ordered K2 fits against exact-three-means projections of their matched
unrestricted fits. It reuses the previous independently qualified witness code
byte for byte, now against the corrected initializer's current results.
It fits nothing and does not open biological or simulated outcome truth.

Only independently verified surviving low-depth full-joint witnesses justify
a subsequent bounded optimizer intervention. Any intervention must pass all
old recovery, 32/64 comparison, false-split/power and uncertainty requirements
before the corresponding promotion. Job **23227019** requests nslab CPU,
1 CPU / 4 GB / 1 h. It is pending with scheduler reason `Priority` at this
checkpoint; no result is assumed.

[Frozen residual-search design](../dev/ordered_residual_search_20261111/DESIGN.md),
[unchanged witness code](../dev/ordered_residual_search_20261111/study.py),
[source identities](../dev/ordered_residual_search_20261111/SNAPSHOT.json),
[submission](../dev/ordered_residual_search_20261111/submission.json).

## Completed residual search and bounded safeguard

The earlier pending statement is historical. Job **23227019** completes 0:0
in 1m13s, scheduler peak **73,908 KiB**. All 5,220 matched comparisons verify.
There are **120** better feasible points, each independently confirmed by full
count PMFs and the same path prior. Among low-depth full joint fits, three
occur with state (indices 205, 207, 216), zero in the 48 clean cases, zero in
48 homogeneous controls and zero in the single state-only control. This is
optimization evidence, not a new population-performance estimate.

The first lexical low-depth witness, 205, has independently reconstructed
objective gain **3.04203278781**, beyond allowance **0.00827556471**. The other
two gains are 25.98347641172 and 0.68124794646, also independently resolved.
The corrected initializer removed most earlier misses but did not eliminate
them. No global optimality certificate follows.

The implemented research safeguard fits the existing unrestricted count
family on each training matrix, projects its matching fit with the existing
exact-three-means rule, and retains a resolved better feasible point under
the unchanged ordered objective. It preserves the original fit and every
optimizer diagnostic, recording when the retained point lacks a new local
convergence certificate. Its cache keys include the exact training matrix;
held-out observations do not initialize those fits. This is an explicit
candidate-set addition, not a claim that more fitting identifies genomic CN.

Preflight **23227108** requests nslab CPU, 1 CPU / 8 GB / 1 h. It runs the
unchanged eight mechanism cases plus witness 205, locks all nine bundles before
truth, independently reconstructs all projected comparisons and invokes the
byte-identical original mechanism reader. Both the original gate and a resolved
same-objective witness improvement are required. The result is pending at this
checkpoint. No full study or public promotion is automatically authorized by
a pass; clean low-depth failures have no corresponding witness in this trace.

[Full residual-search result](../dev/ordered_residual_search_20261111/run_23227019/results/REPORT.json),
[safeguard scientific contract](../dev/projected_incumbent_20261111/DESIGN.md),
[exact safeguard](../dev/projected_incumbent_20261111/fitter.py),
[preflight and independent comparisons](../dev/projected_incumbent_20261111/preflight.py),
[frozen identities](../dev/projected_incumbent_20261111/SNAPSHOT.json),
[exact preflight submission](../dev/projected_incumbent_20261111/submission.json).

## Completed safeguard preflight: better fit, unresolved workflow recovery

Job **23227108** completes 0:0 in 2m26s, with peak process RSS **111,308 KiB**
(scheduler sample 79,480 KiB). All nine bundles lock before outcome truth.
All **162** projected-point comparisons independently verify. Four projected
points are retained, each independently resolved, and their gains propagate
through the existing nested fitting hierarchy. The byte-identical original
reader checks 144 structured mechanism fits and passes the original gate.
All four controls select one population; both high-depth heterogeneous cases
select two; both low-depth heterogeneous cases still select one.

The first witness's full K2 count-MAP objective improves by **6.00835149716**,
beyond numerical allowance **0.00827674258**. Its fixed-K2 state-free minority
MSE falls **0.09407654 → 0.05508351**, but whole state-free MSE rises
**0.00417117 → 0.00432785**. Fixed-K2 whole RNA MSE also rises
**0.00758669 → 0.00817149**, while minority RNA MSE falls
**0.07844567 → 0.04940650**. The report includes conditional quadrature
diagnostics; these are not statistical uncertainty or DNA truth.

The actual joint selector still chooses **K1**, with two-fold predictive score
**-2.8129181e-6** versus numerical diagnostic **0.00237469**. Thus a better
full-sample fitted K2 point is insufficient for a demonstrated workflow recovery
gain. The mechanism pass establishes the bounded computational safeguard and
preserves its controls; it does not repair all seven comparisons. No full
290-case run, new calibration bank or public promotion is submitted.

This result separates two failures: the fitting stage can discard a better
feasible profile, and selection can still withhold that distinction after
fitting improves. The whole/minority risk tradeoff is also real in this witness.
A future diagnostic would have to distinguish reduced training information
from the predictive loss's preference; forcing K2 from the known truth is not
a correction. Original false-split/power and uncertainty criteria remain.

| Job | Completion | Status | Execution errors | Peak memory | Downstream |
|---|---|---|---|---|---|
| 23226980 score/exclusion qualification | 1/1 ████████ | Complete, 0:0 | 0 | 97,028 KiB process | LOO panel not submitted after collapse finding |
| 23227019 residual search | 5,220/5,220 ████████ | Complete, 0:0 | 0 | 73,908 KiB scheduler | Witness justifies bounded safeguard |
| 23227108 safeguard preflight | 9/9 ████████ | Complete, 0:0 | 0 | 111,308 KiB process | No full run submitted |

Running 0, pending 0, completed 3 in this chain. All use account/QOS nslab.
Static checks preserve 91 package/native/data/test identities and all 385/387
prior ledger values at their respective snapshots. Nineteen Python parses,
four shell syntax checks and 836 local links pass at the running-job checkpoint.
The public suite is not rerun because its source remains unchanged. Historical
results and every original scientific criterion remain intact.

[Full independent preflight](../dev/projected_incumbent_20261111/preflight_results/REPORT.json),
[original reader's mechanism results](../dev/projected_incumbent_20261111/preflight_results/mechanism/REPORT.json),
[running-job preservation checkpoint](../dev/projected_incumbent_20261111/preservation_checkpoint_v1.json),
[terminal preservation checkpoint](../dev/projected_incumbent_20261111/preservation_checkpoint_v2.json),
[all three completed scheduler records](../dev/projected_incumbent_20261111/accounting_complete.tsv).

The terminal checkpoint preserves all earlier 385/387 claim values, with 390
entries now recorded, and verifies every before/after execution source log:
1,175 qualification identities, 1,171 residual-search identities and 25 preflight
identities. Both first lexical low-depth full joint contrast vectors are
additionally confirmed exactly zero in the saved parameter records. This
supports the specific collapse diagnosis without broadening it into a
measurement information limit.

Public recovery stays 2/5, the separate terminal-pair research candidate stays
5/5, continuous-state false splits stay 14/200, and selected correlated interval
coverage stays 54/1000. No public default, criterion or frozen outcome changes.
No biological replication or independently qualified genomic truth was added.
