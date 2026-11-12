# Recovery fitting: reuse existing feasible solutions

This follows findings 4 and 12 of the
[scientific review](scientific_code_review_2026-11-10.md). The goal is to repair
a demonstrated optimization inconsistency before evaluating recovery again.
It does not resolve genomic identification or biological error control.

## Cause and fixed intervention

The [completed exact-LOO investigation](LOO_MEMBERSHIP_DIAGNOSTIC_2026-11-11.md)
shows that the final witness score is accurately integrated, while different
fitted solutions produce a large convergence diagnostic. The original failed
criterion is retained. An integration-only diagnostic cannot replace it.

The fixed preflight uses only the training fold excluding cell 91 of witness
205. This was the first lexical maximizer of the original numerical diagnostic,
selected without membership or DNA truth. It first checks whether the existing
96-node fitted point improves the 24-node target over its old point, evaluating
both under that same target. Complete SciPy PMFs and an independently evaluated
path prior must resolve the improvement beyond the existing arithmetic budgets.
Without that witness, no refitting follows.

The candidate passes already evaluated, compatible solutions between the three
original fitting orders. Each training-only sweep preserves its own feasible
point and stops when no native-target gain exceeds the original arithmetic
budget. All steps are reconstructed independently. The 1,000-sweep ceiling
means unresolved if reached. No global-optimum certificate is claimed.

The preflight also requires no independently resolved training-objective loss
at common 384-node precision, and a reduction in the original single-fold joint
diagnostic beyond its numerical error. The excluded cell is scored only after
training stops. Likelihood, prior, dimensions, original optimizer settings and
scientific thresholds remain unchanged. This adds feasible starting points;
it is an optimization change requiring evidence, not a guarantee of recovery.

## Execution and decision

One job, 23227777, is submitted on account/QOS nslab, CPU partition, one CPU,
8 GB RAM and one hour. Static Python parsing, shell syntax and frozen input
identities pass before submission. Its independent same-target witness passes
and training-only refinement starts. The final preflight is pending at this
checkpoint. There are no dependent jobs or arrays.

Only a passing preflight permits a separately frozen rerun of the same three
LOO probes. It replaces none of the original eight-case, 32/64-family,
false-split/power, coverage or biological-validation requirements. Failure stops
expansion. Production recovery remains 2/5, research terminal-pair recovery 5/5,
state false splits 14/200 and selected correlated coverage 54/1,000.

The estimands remain working RNA density and state-zero RNA response. Regional
transcription is still an alternative to genomic dosage. These are exposed
simulations, biological n=0; independent corresponding DNA and suitable
membership-resolution truth remain necessary for stronger claims. No production
code, frozen output or historical claim is changed by this preflight.

[Frozen design](../dev/shared_incumbent_20261111/DESIGN.md),
[implementation](../dev/shared_incumbent_20261111/model.py),
[independent checks](../dev/shared_incumbent_20261111/preflight.py),
[source snapshot](../dev/shared_incumbent_20261111/SNAPSHOT.json),
[submission](../dev/shared_incumbent_20261111/submission.json).

## Completed one-fold result

Job 23227777 completes 0:0 in 23 scheduler seconds, 22.72 process seconds,
137,556 KiB process peak RSS. Both source checks pass. The independent witness
gain is 3.638736926 against budget 0.001235187. Refinement stabilizes after two
sweeps; all 24 independently reconstructed steps and all 12 common-384-node
nondegradation comparisons pass. Joint single-fold diagnostic falls from
4.277144984 to 0.000150272, exceeding the numerical bound on that reduction.
No-factor diagnostic is unchanged. Fit-96 held-out score remains 2.400642842;
the correction primarily repairs earlier-order solutions in this fold.

This passes the declared preflight. It does not establish that the summed
120-fold diagnostic resolves, that controls retain one component, or that
actual profiles improve. Only the separately frozen same-three-probe test is
justified next. Earlier pending statements above are historical checkpoints.

[Independent result](../dev/shared_incumbent_20261111/run_23227777/results/REPORT.json),
[all update records](../dev/shared_incumbent_20261111/run_23227777/results/STABILITY.json),
[process resources](../dev/shared_incumbent_20261111/run_23227777/resources.txt),
[scheduler accounting](../dev/shared_incumbent_20261111/accounting_preflight.tsv).

## Frozen same-three-probe continuation

The passing preflight permits exactly the declared follow-up, now isolated in
`dev/shared_incumbent_loo_20261111`. Writers 23227839, 23227840 and 23227841
refine the original 360 training-only cached families. Reader 23227842 depends
on all three succeeding. All use account/QOS nslab, CPU, one CPU; writers have
8 GB and two hours, the reader 4 GB and two hours. No arrays are submitted.

All original reader function bodies and the scoring/selection function ASTs
are exact. Only the fitting function changes, with per-fold IDs and training
hashes checked against the locked old points. Each case repeats the excluded-
only perturbation test. All update histories are saved and independently
reconstructed before the reader calls its original target/scoring functions.
Every new output locks before any truth construction. The same old descriptive
full-data fits isolate the effect of selection. There are 62 frozen source/input
identities checked before and after execution. Outcomes are pending here.

[Frozen follow-up design](../dev/shared_incumbent_loo_20261111/DESIGN.md),
[fitting implementation](../dev/shared_incumbent_loo_20261111/selection.py),
[additional independent checks](../dev/shared_incumbent_loo_20261111/verify_updates.py),
[source snapshot](../dev/shared_incumbent_loo_20261111/SNAPSHOT.json),
[writer submissions](../dev/shared_incumbent_loo_20261111/case_jobs.json),
[reader and dependencies](../dev/shared_incumbent_loo_20261111/submission_reader.json).

## Completed three-probe result: numerical repair restores the distinction

All four jobs complete 0:0 on nslab. All 360 excluded-cell families stabilize;
5,184 independently reconstructed native-target updates and 4,320 common-384
old/new comparisons pass. All 2,160 full-PMF pairs, 720 higher-order pairs,
excluded-only perturbations and exact unchanged decisions pass. Every output
locks before target construction. All 62 before/after identities pass in each
of the four executions. There is no restart or changed criterion.

| Joint route | Selected components | Final score | Original numerical diagnostic |
|---|---:|---:|---:|
| Homogeneous control | 1 | 0.000002888 | 0.000091255 |
| Recovery witness 205 | 2 | 6.324223511 | 0.008844259 |
| State-only control | 1 | -0.000007865 | 0.000121332 |

The original three-probe mechanism criterion passes. Its preceding LOO score
5.852086 was below diagnostic 36.786157; that failed result remains unchanged.
The correction reduces fitted-solution inconsistency under the same targets
and makes the witness decision resolvable. No integration-only diagnostic is
substituted. The no-factor route still splits the state-only control and is
not an acceptable genomic rescue.

Witness globally matched MAP assignments improve from 96/120 to **105/120**.
The smaller fitted component has **17.880703 responsibility mass and 17 MAP
cells**, below the unchanged public minimum 30. This is a working RNA mixture
distinction, **not a supported public clone call**. No small-component criterion
is relaxed. Both controls remain one group with 120/120 matched assignments.

| Witness point MSE | Prior actual K1 | New selected K2 | Pseudobulk | Empirical |
|---|---:|---:|---:|---:|
| Whole RNA | 0.00758669 | 0.00817149 | 0.01448416 | 0.01113565 |
| Minority RNA | 0.07844567 | 0.04940650 | 0.12567801 | 0.06558275 |
| Whole state-zero response | 0.00417117 | 0.00432785 | 0.00442919 | 0.01025733 |
| Minority state-zero response | 0.09407655 | 0.05508351 | 0.08742357 | 0.06835610 |

The minority gain and adverse whole-profile point change are both retained.
The unchanged descriptive full-data fits still differ across fitting orders;
their conservative whole-response numerical allowance is 0.00261211 and
minority allowance 0.01061304. These are numerical diagnostics, not confidence
intervals. The small whole-response difference from pseudobulk is unresolved
under that allowance. Do not call the point means a broad accuracy pass.

The next justified correction is to apply the same qualified optimizer to the
three full-data descriptive families, after preserving independent selection
statistics. It may resolve numerical ambiguity in profiles; it cannot change
the held-out decision or waive the component support minimum. No original
eight-case, 32/64-family, production 2/5, research terminal-pair 5/5, false-split,
power, uncertainty or biological-validation outcome is replaced. Biological
n remains zero and no SOTA claim follows.

| Job | Work complete | State / exit | Elapsed | Scheduler peak RSS |
|---|---|---|---|---:|
| 23227839 | 120/120 ████████ | Complete, 0:0 | 17m45s | 66,036 KiB |
| 23227840 | 120/120 ████████ | Complete, 0:0 | 19m08s | 66,336 KiB |
| 23227841 | 120/120 ████████ | Complete, 0:0 | 12m21s | 64,328 KiB |
| 23227842 | 360/360 ████████ | Complete, 0:0 | 4m23s | 134,700 KiB |

Running 0, pending 0, done 4, execution errors 0 in this chain; no downstream
job is submitted at this checkpoint. Earlier pending statements are historical.

[Complete independent result](../dev/shared_incumbent_loo_20261111/reader_23227842/results/REPORT.json),
[all independent training checks](../dev/shared_incumbent_loo_20261111/reader_23227842/results/UPDATE_CHECKS.json),
[complete accounting](../dev/shared_incumbent_loo_20261111/accounting_complete.tsv).

## Same correction for descriptive full fits

The three-probe pass justifies applying the identical qualified correction to
the three descriptive full-data families. Job **23228204** uses nslab CPU,
1 CPU, 8 GB and one hour. Its final outcome is pending at this
checkpoint. Independent selection statistics and all previous results remain
unchanged. All native objective updates, common-384 targets and posterior
profiles are checked independently; all RNA outputs lock before truth scoring.
The original profile-order diagnostic remains, with 384 nodes only checking
integration. Small-component support below 30 still withholds public clone
support. No original scientific gate or false-split result is replaced.

[Frozen refit design](../dev/shared_incumbent_refit_20261111/DESIGN.md),
[executed implementation](../dev/shared_incumbent_refit_20261111/study.py),
[source snapshot](../dev/shared_incumbent_refit_20261111/SNAPSHOT.json),
[exact submission](../dev/shared_incumbent_refit_20261111/submission.json).

## Completed full-fit result: profile comparisons now resolve numerically

Job23228204 completes0:0 in43 scheduler seconds. The three families stabilize
after1/2/1 sweeps. All48 independently reconstructed objective updates,36
common384 nondegradation checks and180 native/independent posterior-array
comparisons pass. All74 before/after source identities pass. Every held-out
selection is exact and all three new RNA outputs lock before target scoring.

The witness's original maximum state-zero profile-order diagnostic falls from
0.15021955 to0.0000974163; its RNA counterpart falls0.10013271 to0.0002435711.
Fixed-fit192-to384 discrepancies are at most3.775e-14. The old diagnostic
definition remains unchanged. The full-data fit96 point barely changes;
the correction chiefly brings its earlier-order fits into agreement. This
does not improve a reported statistical confidence interval.

| Witness endpoint | Stabilized point MSE | Numerical allowance | Pseudobulk | Prior actual K1 |
|---|---:|---:|---:|---:|
| Whole RNA |0.0081714893|0.0000096983|0.0144841639|0.0075866927|
| Minority RNA |0.0494065043|0.0000282135|0.1256780073|0.0784456701|
| Whole state-zero response |0.0043278548|0.0000058624|0.0044291862|0.0041711740|
| Minority state-zero response |0.0550835059|0.0000289111|0.0874235693|0.0940765547|

The pseudobulk differences now exceed their numerical allowances. This is a
specific exposed simulation gain, not an independent biological superiority
result or a replacement for the failed full32/64-family gate. The adverse
whole-profile change versus priorK1 also exceeds numerical uncertainty and
remains a real tradeoff. Every previous empirical/other baseline stays in the
report, including their complete earlier endpoint records.

Globally matched assignments remain105/120 after the full refit. Smaller
component mass is17.880703 with17 MAP cells, still below public30. No supported
genomic clone is emitted. Both controls retainK1 and120/120 matched assignments;
homogeneous joint RNA MSE remains worse than pseudobulk. No-factor still splits
the state-only control. These adverse observations prevent an all-task benefit
claim. Biological n=0, genomic attribution and biological error control remain
unestablished. Production recovery2/5,state errors14/200 and selected coverage
54/1000 remain unresolved and unchanged.

[Full numerical and performance report](../dev/shared_incumbent_refit_20261111/run_23228204/results/REPORT.json),
[RNA lock before targets](../dev/shared_incumbent_refit_20261111/run_23228204/results/RNA_LOCKED.json),
[scheduler accounting](../dev/shared_incumbent_refit_20261111/accounting_complete.tsv),
[process resource record](../dev/shared_incumbent_refit_20261111/run_23228204/resources.txt).

## Preservation and remaining work

The final checkpoint verifies91 protected package/native/data/test identities,
all394/396/398 prior ledger values at their respective task baselines, seven
Python sources, four shell launchers,868 existing local report links and every
before/after execution-source log. The ledger has400 append-only entries.
The preflight and all five continuation jobs complete0:0 on nslab CPU; there
are no remaining task jobs or downstream dependencies. The public suite is
not rerun because its validated source identities remain exact.

This is a real numerical/selection correction with a specific minority-profile
and membership gain on one exposed probe. It does not fix all collapsed low-depth
fits, the three original production recovery failures, complete-workflow
false-split/power, selected uncertainty or missing independent genomic truth.
Those original requirements remain. Further fitting changes need a demonstrated
objective or information-loss mechanism; no arbitrary start/rank/threshold sweep
or consumed-panel retuning is justified. The full genomic/clone/spatial/SOTA goal
remains active and unachieved.

[Final preservation checkpoint](../dev/shared_incumbent_refit_20261111/preservation_checkpoint_v1.json).
