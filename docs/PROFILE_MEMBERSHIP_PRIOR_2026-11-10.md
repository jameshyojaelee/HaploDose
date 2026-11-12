# Discovery profile deviations and clone prevalence

The new regional RNA profile has not yet demonstrated better clone recovery.
This experiment tests whether its discovery-only cell deviations improve the
existing conditional comparison between two RNA populations and continuous
expression state. The implementation and nine input partitions are locked.
The initial replay failures are resolved by matching each case's archived
qualified runtime: all 27 arrays and scores now match exactly. The full screen
rejects both candidates, at 6/9 for balanced prevalence and 7/9 for the profile
partition. Balanced prevalence restores positive conditional evidence for all
five original contrasts but also accepts a state false distinction. Production
recovery remains 2/5; no full-workflow correction or error guarantee follows.

## Question and prespecified decision

A common regional gene-profile anchor is constant across cells. It cannot
itself improve their separation within that region. The proposed initializer
therefore transfers only cross-fitted expected-log RNA fractions relative to
the archived expected fractions. It does not transfer a DNA label or an
uncalibrated probability as a genomic membership probability.

The panel remains all five original recovery cases and four previously
identified expression-state false distinctions. For each saved parent and
selected interval, aggregate the proposal gene counts into the original bins;
require exact identity alignment and count conservation. Divide cells into
three seeded, identity-aligned folds. Estimate variable-library beta moments
from the other folds, calculate expected log fractions, subtract the archived
expected log fractions, apply the fixed five-bin median within the interval,
and average. An exact rational two-mean partition determines discovery
prevalence in each existing depth stratum. Constant scores, nonunique optima,
small groups or absent within-stratum component support abstain.

Two fixed arms compare these prevalences with the simple control of one-half
in every stratum. Both retain the archived gene NB observation model, expected
counts, sizes, depth strata, optimizer and 24/48/96-point quadrature. New
profiles are locked before the historical confirmation replay; new fits must
run separately and freeze both models before confirmation counts are opened.
The historical nuisance estimation and interval selection remain conditioned
upon. This is not a complete independent-confirmation workflow.

Every original contrast must have positive resolved predictive advantage,
and every state contrast negative resolved advantage: 9/9 per arm. The old
numerical-resolution diagnostic is unchanged. Unsupported cases fail their
arm. All 27 historical per-cell prediction arrays and final scores must
first replay exactly. Only a complete pass permits a separately frozen
full-workflow recovery and fresh false-split/power study. The 5% criterion,
confidence bounds and power assessment do not change. If both arms pass,
prefer the simpler balanced arm absent separate evidence for the profile arm.

## Executed checks and failed attempt

Job 23191185 used nslab account/QOS, CPU partition, one CPU, 8 GB and a
one-hour limit on node003. Installed package/native identities and the
source checks passed. The following independent checks passed before replay:

- Variable-depth moment identities and exact target-fold isolation.
- The known Beta(2,3) expected-log value and constant-score abstention.
- Exact exhaustive binary-allocation agreement for the two-mean fixture,
  identity-aligned permutation, and unchanged-peer exclusion.
- Existing rational NB probabilities, independent full-PMF quadrature,
  heat-equation gradients, nested-null and component/stratum permutations.

All nine proposal profiles were saved before new confirmation scoring;
all nine satisfy the discovery partition support requirements. This is
input qualification, not nine successful recovery results.

The preflight then stopped at exact historical-array equality: 22/600
entries differ, with maximum absolute discrepancy 4.440892098500626e-16.
No parity completion file was written, so candidate fitting remains
ineligible. The failed loop does not record the case index. [Inference]
The discrepancy is consistent with the earlier node003 case-06 replay
failure, whose unchanged retry passed on node010; this does not by itself
prove the cause of the current failure.

Job 23191239 repeats the same source, environment and exact criterion on
node010, writing a fresh `prepared_v2` destination. It is pending with
reason `Resources` at this checkpoint. No running job was disturbed.
The first submission-shell attempt used unavailable `python` and exited
before submitting; using `python3` submitted the retry once. There is no
candidate fit job or downstream summary submission.

| Job | Work | Progress | State | Memory | Downstream |
|---|---|---|---|---|---|
| 23191185 | Input/numerical checks and exact replay | 0/1 ░░░░░░░░ | Failed, 1:0, 34 s | Scheduler RSS 120,256 KiB; process peak 185,744 KiB | Candidate fits blocked by failed exact gate |
| 23191239 | Unchanged node010 replay | 0/1 ░░░░░░░░ | Pending: Resources | 8 GB requested; usage not available | Fits not submitted |

Current batch: zero running, one pending, zero successfully completed,
one failed. Numerical checks passing within a failed job do not make that
job successful. The failed outputs, source hashes and previous evidence
are retained. No production code, test, frozen artifact or criterion changed.

## Scientific interpretation and evidence still needed

The estimand is conditional RNA predictive advantage of a two-population
mixture over continuous expression state. Working beta moments initialize
prevalence; shared-state gene NB scores assess it. These stages are not an
exact joint model. Regional transcription, reference mismatch and unmodelled
nuisance can still reproduce a CNA-shaped signal. Biological n is zero;
selected simulations supply no biological false-split guarantee.

A failed nine-case screen rejects this initializer without threshold tuning.
A passing screen would support only workflow integration. Full-workflow
null/power experiments, corresponding-cell DNA for membership, and independent
single-CNA biological specimens at documented resolution remain necessary.
Neither more quadrature nor an exact replay supplies missing replication.

[Prespecified design](../dev/profile_membership_prior_20261110/DESIGN.md),
[implementation](../dev/profile_membership_prior_20261110/study.py),
[profile calculation](../dev/profile_membership_prior_20261110/profiles.py),
[source hashes](../dev/profile_membership_prior_20261110/source.sha256),
[locked inputs and numerical checks](../dev/profile_membership_prior_20261110/prepared_v1/DISCOVERY_INPUTS_LOCKED.json),
[failed replay](../dev/profile_membership_prior_20261110/analysis_23191185/output.txt),
[original submission](../dev/profile_membership_prior_20261110/prepare_submission.json),
[unchanged retry](../dev/profile_membership_prior_20261110/prepare_retry_submission.json),
[scheduler checkpoint](../dev/profile_membership_prior_20261110/scheduler_checkpoint.txt),
[accounting checkpoint](../dev/profile_membership_prior_20261110/accounting_checkpoint.tsv),
[earlier platform replay](FIXED_DOSAGE_STATE_2026-11-10.md).

## Resolved replay: numerical environment, not changed observations

The initial node010 retry 23191239 fails in the same way as node003. Diagnostic
23194965 then verifies identical input values and repeats all 27 predictions
with original, serialized, C-order and F-order arrays. Only case06's three
arrays differ in every layout. At orders 24/48/96 the mismatch counts are
22/24/22, with maximum discrepancies 4.44e-16/6.66e-16/4.44e-16. Memory order
and serialized input values do not explain these differences.

The earlier successful case06 retry used `replay_native.sbatch`, not the copied
runtime launcher. It invokes the existing Conda interpreter with environment
libraries first, allowing the node's native libc/libm. The copied runtime003
loader instead matches the eight other cases. A native attempt on node003,
23195230, fails case06. It inherited node003 from its header despite an announced
node010 intention; it was left running and retained. Native node010 attempt
23195367 instead fails already at case00. The actual mapped library paths/hashes
are recorded; the two nodes have different libc, libm and loader identities.
The particular libm instruction responsible has not been isolated.

Every decompressed NPY entry in all nine prepared inputs is byte-identical
across copied and native-node010 executions, including each candidate partition.
The archived fixed-dosage study passed with copied runtime003 for eight cases
and native node010 for case06. A single-process preflight using either of those
environments need not reproduce that whole panel. This is a runtime-reproduction
problem, not a changed biological observation or a failed recovery outcome.

Qualification 23195868 repeats case06 in native node010 and verifies the other
24 exact comparisons from the full copied-runtime diagnostic. It requires every
original array and every order-specific/final score, verifies byte-identical
new inputs and writes a fresh `prepared_v5`. **All 27 exact comparisons pass.**
The equality, observations, candidate partitions, optimizer, numerical diagnostics
and scientific decisions are unchanged. Earlier failures remain failures.
Later workers use their case's qualified runtime and repeat their exact checks.
The preflight now supports multiple recorded environments; no tolerance or
case omission replaces any archived comparison.

| Job | Attempt | Progress | Outcome | Requested resources |
|---|---|---|---|---|
| 23191185 | Copied runtime, node003 | 0/1 ░░░░░░░░ | Failed exact check, 34 s | nslab CPU1 / 8 GB / 1 h |
| 23191239 | Copied runtime, node010 | 0/1 ░░░░░░░░ | Failed exact check, 36 s | nslab CPU1 / 8 GB / 1 h |
| 23194965 | Values and four memory layouts | 9/9 ████████ | Completed diagnostic, 30 s; did not pass case06 | nslab CPU1 / 4 GB / 15 min |
| 23195230 | Native interpreter, inherited node003 | 0/1 ░░░░░░░░ | Failed exact check, 16 s | nslab CPU1 / 8 GB / 1 h |
| 23195367 | Native interpreter, explicit node010 | 0/1 ░░░░░░░░ | Failed exact check at case00, 7 s | nslab CPU1 / 8 GB / 1 h |
| 23195868 | Case-specific runtime qualification | 27/27 ████████ | Completed, 16 s | nslab CPU1 / 4 GB / 15 min |

The nine case jobs 23195991–23195999 each request nslab CPU1 / 8 GB / 2 h,
without an array. Case06 uses native node010; the other cases use copied
runtime003 on node010. Summary 23196000 requests CPU1 / 4 GB / 15 min and
depends on all nine. At this checkpoint eight cases complete, one is running,
and the summary is pending. No candidate result is promoted from partial cases.
Package/native source and existing scientific tests remain unchanged.

[All input files byte-identical](../dev/profile_membership_prior_20261110/prepared_runtime_file_comparison.json),
[full memory/runtime diagnostic](../dev/profile_membership_prior_20261110/diagnosis_23194965/results/REPORT.json),
[native node003 library identities](../dev/profile_membership_prior_20261110/analysis_23195230/identity.json),
[native node010 library identities](../dev/profile_membership_prior_20261110/analysis_23195367/identity.json),
[case-specific qualification contract](../dev/profile_membership_prior_20261110/CASE_RUNTIME_QUALIFICATION.md),
[qualification implementation](../dev/profile_membership_prior_20261110/qualify_replay.py),
[27 exact matches](../dev/profile_membership_prior_20261110/prepared_v5/BASELINE_PARITY.json),
[fresh locked inputs](../dev/profile_membership_prior_20261110/prepared_v5/DISCOVERY_INPUTS_LOCKED.json),
[all case submission headers and commands](../dev/profile_membership_prior_20261110/case_submissions.jsonl),
[dependent summary](../dev/profile_membership_prior_20261110/summary_submission.json).

## Completed nine-case screen and decisive tradeoff

All nine workers and summary 23196000 complete `0:0`. Each worker reproduces
its three original per-cell arrays and final score before recording a new
scientific result: all 27 exact checks pass again. No arm or case disappears.
Models for both arms are locked before new confirmation scoring.

The following are conditional mixture-versus-continuous RNA predictive
advantages. Positive resolved values favour the mixture; negative resolved
values favour continuous state. `Unresolved` fails the locked criterion and
is not reported as a false positive or a successful rejection.

| Locked case | Original soft prevalence | Balanced prevalence | Profile partition prevalence |
|---|---:|---:|---:|
| Original 20261031 | 4.348085 | 4.684102 | 2.245880 |
| Original 20261107 | 4.713200 | 0.019035 | 1.295589 |
| Original 20261108 | 4.527263 | 4.226870 | 0.310007 |
| Original 20261109 | -3.396705 | 0.497250 | -1.16e-7, unresolved |
| Original 20261110 | 2.863255 | 0.891984 | 1.425332 |
| State 202611096001 | -2.277543 | 5.11e-7, unresolved | -2.137681 |
| State 202611096002 | -1.836996 | 2.45e-7, unresolved | -0.620129 |
| State 202611096004 | -4.167723 | 0.605477, false distinction | 2.426788, false distinction |
| State 202611096006 | -10.889296 | -1.115520 | -4.491381 |

The original native-q screen remains 8/9. Balanced prevalence passes all five
original contrasts but only one state control, with two unresolved controls
and one resolved false distinction: 6/9 overall. The new profile partition
passes four original contrasts and three state controls: 7/9 overall.
The difficult original contrast remains unresolved, and the same state case
is accepted. **Reject both new arms at the unchanged 9/9 continuation gate.**
Do not integrate either into the full workflow or expand its null/power study.

This changes the mechanism assessment. Holding observations and likelihood
fixed, prevalence assumptions alter the difficult recovery's conditional
support. They do not add genomic information. Balanced prevalence restores
that contrast while losing state specificity; the RNA profile partition does
not repair it. Membership initialization alone does not resolve this tradeoff.
These results do not prove that RNA from the cases is intrinsically incapable
of recovery under every alternative method.

Threshold tuning cannot repair this statistic on the fixed panel. Under the
balanced arm, the false state score 0.605477 exceeds the true 20261109 score
0.497250 and true 20261107 score 0.019035. Under the profile arm, the false
score 2.426788 exceeds even the largest recovered original score, 2.245880.
Any higher scalar cutoff that removes those false distinctions also removes
required true contrasts. The ordering exceeds the existing numerical
diagnostics; those diagnostics remain convergence checks, not formal
integration enclosures or population error intervals.

Independent reader 23198336 imports none of the fitting/scoring code. It verifies
the exact nine identities and order, computes all **54** saved per-cell sums
using exact rational arithmetic, independently recounts all **18** numerical
resolution/criterion decisions, and checks the cutoff-obstruction witnesses.
All checks pass. This is computational verification of the recorded negative
decision, not biological validation. These selected simulation cases have
biological n=0 and do not estimate a specimen false-split rate.

| Job | Case / work | Progress | Status and elapsed |
|---|---|---|---|
| 23195991 | Original 20261031 | 1/1 ████████ | Completed, 3m21s |
| 23195992 | Original 20261107 | 1/1 ████████ | Completed, 1m59s |
| 23195993 | Original 20261108 | 1/1 ████████ | Completed, 2m13s |
| 23195994 | Original 20261109 | 1/1 ████████ | Completed, 3m24s |
| 23195995 | Original 20261110 | 1/1 ████████ | Completed, 3m53s |
| 23195996 | State 202611096001 | 1/1 ████████ | Completed, 4m30s |
| 23195997 | State 202611096002, native runtime | 1/1 ████████ | Completed, 6m09s |
| 23195998 | State 202611096004 | 1/1 ████████ | Completed, 2m49s |
| 23195999 | State 202611096006 | 1/1 ████████ | Completed, 2m21s |
| 23196000 | Complete fixed-family summary | 9/9 ████████ | Completed, 13s |
| 23198336 | Independent rational recount | 54/54 ████████ | Completed, 4s |

The recount uses nslab CPU1 / 4 GB / 15 min on v4 Intel; case and summary
resources are above. Across all 17 experiment jobs: 13 complete, four failed
preflights, zero running and zero pending. All four failures remain. No full
package/native/GPU suite was rerun for these research-only changes. No
production code, frozen result or scientific threshold changed. The original
2/5 production recovery, 14/200 continuous-state false-split result, 195/200
power result and selected-interval uncertainty failure remain separate evidence.

[All final scores, diagnostics and decisions](../dev/profile_membership_prior_20261110/summary_v1/report.json),
[independent exact recount](../dev/profile_membership_prior_20261110/recount_23198336/results/REPORT.json),
[independent implementation](../dev/profile_membership_prior_20261110/recount.py),
[recount source identities](../dev/profile_membership_prior_20261110/recount_source.sha256),
[complete submission journal](../dev/profile_membership_prior_20261110/jobs_complete.json),
[final scheduler accounting](../dev/profile_membership_prior_20261110/accounting_complete.tsv),
[all measured process resources](../dev/profile_membership_prior_20261110/resource_summary_complete.json),
[125 manifest entries verified](../dev/profile_membership_prior_20261110/source_integrity_complete.json),
[prior ledger values preserved](../dev/profile_membership_prior_20261110/ledger_complete_preservation.json),
[original native-q comparison](DEPTH_STATE_PREDICTIVE_2026-11-10.md).
