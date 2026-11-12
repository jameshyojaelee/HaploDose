# Is the remaining failure in the alternative fit or its selection?

The [completed excluded-cell experiment](INTERVAL_INCUMBENT_LOO_2026-11-11.md)
improves one clean simulated dosage-component profile, but leaves the
state-confounded case at K1. Its complete criterion fails. This bounded trace
uses the already frozen selected and unselected models to locate the remaining
loss before proposing any correction. No new estimator is fitted.

## Fixed measurements and decisive observations

Keep all four original cases 0, 2, 4 and 6, their full K1/K2 models, both state
branches, parameters, identities, counts and recorded selection. Verify every
source prediction lock before reconstructing the same simulated targets.
Evaluate all full models at common 384-node precision with independent SciPy
full PMFs, 192/384 convergence and adaptive probes. Keep original 24/48/96
fitted-point sensitivity when reporting numerical endpoint allowances.

Report whole/minority RNA and state-zero errors, matched membership, component
support, paths and training-MAP contrasts for every model. Decompose the existing
held-out score into simulated minority and majority contributions, with exact
Fraction conservation of the total. Centred loading/contrast alignment to the
known simulated dosage direction is descriptive geometry, not causal proof.

If frozen K2 improves both dosage-component endpoints beyond numerical
uncertainty despite negative held-out prediction, the selection score does not
preserve that target in this case. If K2 is also inaccurate, accepting it cannot
establish a repair. Mixed or numerically unresolved results stay explicit.
Neither result establishes global optimality or an RNA information limit.
Truth never chooses a model, parameter, region or threshold for deployment.

The estimands remain conditional RNA prediction and the known simulated
state-zero RNA response. Coordinated transcription can mimic genomic dosage;
independent corresponding DNA and suitable membership truth remain necessary.
Biological n=0. No public calls, original criterion, consumed-dataset role,
frozen result or historical decision changes. The failed full mechanism and
the original recovery/error/power/coverage requirements remain binding.

## Execution record

Job **23230235** requests nslab/nslab CPU, 1 CPU, 4 GB and 1 hour. It uses the
existing isolated read-only Conda/runtime, one thread, neutral working directory
and 247 before/after identities. No arrays, new dependencies, shared installs,
public edits or unrelated-job changes. Static Python/shell checks pass, including
an AST check that the new script calls no fitting or minimization method.
An unmatched parenthesis was caught and corrected before freezing or submission;
no failed Slurm attempt resulted. The completed outcome follows.

[Frozen design](../dev/frozen_population_trace_20261111/DESIGN.md),
[complete diagnostic](../dev/frozen_population_trace_20261111/study.py),
[source identities](../dev/frozen_population_trace_20261111/SNAPSHOT.json),
[exact submission](../dev/frozen_population_trace_20261111/submission.json).

## Completed result: accepting the alternative would not repair this failure

All 16 frozen models were evaluated without refitting. Independent full PMFs,
192/384-node convergence, adaptive probes and exact Fraction score conservation
pass. All 247 before/after source identities and all 91 protected package,
test and native-source identities agree. This is a computational diagnostic,
not a fresh execution of the public test suite or a biological validation.

For the state-confounded case4, retaining the continuous state gives:

| State-zero RNA endpoint | Selected K1 | Rejected K2 |
|---|---:|---:|
| Whole-profile MSE | 0.00386250 | 0.00449645 |
| Minority-profile MSE | 0.08189761 | 0.07913574 |
| Globally matched cells, all 120 included | 96/120 | 93/120 |
| Smaller responsibility mass | No second component | 8.83281 |

The whole-error increase exceeds the summed numerical allowances
0.00001223 + 0.00002301. The minority improvement also exceeds its numerical
allowances, but still loses to pseudobulk's 0.07513644. These allowances retain
original fitted-point sensitivity; they are not statistical confidence bounds.
The alternative's smaller mass is below the unchanged public minimum30.

Both simulated populations lose excluded-cell predictive score: majority96
contributes -1.799716 and minority24 contributes -2.602052, exactly conserving
the existing total -4.401768. Thus this rejection is not simply majority loss
outweighing positive minority predictive support. The full alternative also
has worse penalized training fit: advantage -3.620677, arithmetic allowance
0.045555. Its broad half-block path and modest alignment to the generating
contrast describe this fitted point; they do not establish the cause of the
miss, global optimality, or an information limit.

Removing the continuous state is contradicted by the controls. In case6,
the no-state K2 receives score +111.881918 for a single simulated CNA population
with expression state. The state-enabled K2 has zero contrast and no resolved
second profile. In case4, the no-state K2 matches only64/120 cells. Removing
state adjustment or forcing a two-group answer is therefore not a justified
performance correction.

The clean case2 retains its partial gain: selected K2 state-zero whole/minority
MSE 0.00287302/0.06558946 and104/120 matches. Majority score -1.115953 plus
minority +1.834945 exactly conserves the positive total +0.718992. Its mass
9.50550 remains below30. Whole RNA prediction still loses to pseudobulk.
Both state-enabled controls remain selected K1. The complete prior mechanism
criterion still FAILS; no scientific criterion or genomic claim changes.

**Decision.** Reject a selector-only rescue and removal of the state component.
The next fitting intervention needs an independently verified better feasible
point or a specific model failure. Do not expand the case bank, retune the
prior, lower support30, or reuse oracle parameters as deployable information.
The existing [actual-count known-law study](JOINT_INFORMATION_ORACLE_2026-11-11.md)
already establishes that parameter knowledge can improve these simulated
profiles; another latent-rate oracle replication would not diagnose this fit.

Job23230235 completes0:0 in1m07s; timed process1m06.25s, process peak106688KiB,
scheduler batch peak66644KiB. The difference reflects scheduler sampling.
Running0, pending0, done1, failed0 for this experiment; no downstream job.
[Authoritative result](../dev/frozen_population_trace_20261111/run_23230235/results/REPORT.json),
[process resources](../dev/frozen_population_trace_20261111/run_23230235/resources.txt).
All pre-existing ledger values and frozen artifacts remain unchanged. Original
public recovery2/5, state false splits14/200, power195/200 and selected correlated
coverage54/1000 remain unresolved. Biological n=0; the full goal is unachieved.
