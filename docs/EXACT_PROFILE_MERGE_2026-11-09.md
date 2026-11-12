# Exact RNA-profile redundancy repaired without a similarity threshold

Expression-only, nonspatial fits now coalesce exactly equal finite tumour dose
vectors before hard-MAP eligibility can hide them from pairwise merging. The
normal component is preserved. Allele-enabled and spatial fits keep their
existing behavior. This is a computational correction, not a claim that genomic
clones with similar RNA profiles are biologically equivalent.

For identical component expression likelihoods f, replacing pi_a*f+pi_b*f with
(pi_a+pi_b)*f leaves the observed mixture likelihood unchanged. Responsibilities
and pooled counts add. The implementation uses exact array equality, no closeness
tolerance or new scientific cutoff. A purely algebraic collapse retains its
fitted rates; a subsequent statistical merge still refits pooled profiles.
The existing search then continues and performs its terminal merge checks.

The identifying assumption is identical expression likelihoods under the shared
reference, nuisance and programme parameters. Additional allele or spatial terms
can distinguish equal expression profiles, so the correction is disabled in those
modes. Independently measured DNA could reveal distinct clones with identical
RNA effects; RNA alone would not resolve them. The retained output is an RNA
profile, not certified genomic identity or calibrated membership confidence.

## Executable evidence

The initial dirty source is preserved under
`dev/exact_profile_merge_20261109/baseline`. Build/test source is `source_v2`.
The noneditable wheel is built and installed in the task's private Conda prefix,
then tested from a neutral directory. Job23121461 passes403 Python tests and
four Rust tests; one GPU-only test is skipped. The wheel SHA256 is
`2fbd5df4ffa76505e740bf4b3e2de58e913051a3d84265860e74c0925f3c14ff`.
The native SHA256 remains
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.

Eight new fixtures check a four-identical-component counterexample, marginal
likelihood and posterior/count conservation, normal identity, near/unequal and
nonfinite profiles, allele/spatial isolation, and exact collapse followed by an
ordinary statistical merge. Floating-point allowances scale with operation
counts and machine precision, not scientific decision thresholds. One existing
synthetic model fixture receives its missing `alleles=None` attribute; its
assertions and all original recovery criteria are unchanged.

The six locked mechanism replays and independent audit pass in
`mechanisms_v3` and `mechanisms_audit_v3` (jobs23121721/23121722). Instrumented
tables and posterior arrays exactly match the first corrected replay. Each
captured algebraic merge is independently reconstructed from exact dose-vector
groups; mixture weights, aggregate responsibilities and normalized rate mixtures
agree within operation-based budgets. No final exact tumour duplicates remain.

| Condition / replicate | Frozen tumour profiles | Corrected profiles | Correct assignments /600, before → after |
|---|---:|---:|---:|
| Matched /0 |1|1|600 →600|
| State /0 |1|1|599 →599|
| State /2 |2|2|499 →499|
| State /1 |4|1|600 →600|
| Two-clone /0 |2|2|453 →453|
| Two-clone /47 |3|2|441 →445|

Four cases have all five public tables exactly unchanged. The state /1 search
now removes redundant profiles earlier and never constructs the old four-profile
terminal configuration. Two removals occur along the revised trajectory, not
three removals at the old terminal state. Two-clone /47 changes after an earlier
exact redundancy is removed; its old final near-equal profiles are not directly
collapsed by an approximate rule. State /2's unequal state-driven distinction
remains, so this correction cannot establish the5% scientific target.

Original production recovery remains2/5 (three failures unchanged); the
effective-unit research path retains5/5 with its known gene-size oracle. Both
are independently recounted under all30 original criteria, and the research
veto receives scalar/exact-sign auditing. Jobs23121565/23121602 and
23121566/23121603 pass. No deployable resolution of the original three failures
is claimed.

## Failed audits retained and corrected

Audit23121601 incorrectly expected three duplicate removals, ignoring how earlier
merging changes the search trajectory. The replacement checks actual algebraic
invariants instead. Audit23121687 then found incomplete profiling: the event
proposal's temporary profiler cleared the enclosing merge profiler. The research
adapter now restores the previous profiler. Replay23121721 establishes exact
output parity, and audit23121722 passes every captured merge invariant. Both
failed audit outputs remain preserved. No production or scientific test is
relaxed to accommodate these instrumentation errors.

The superseded pending build23121367 was cancelled before execution; staging
23121362 and replacement full build23121461 pass. All commands, environments,
resources and dependency records are in `jobs.json`, `analysis_JOB/` and
`validation_JOB/` within the task directory.

## False-split, power and uncertainty gate

`PANEL_DESIGN.md` fixes a paired600-fit replay of the previously exposed panel:
100 generating seeds across three conditions and raw/effective research units.
This is regression evidence, not independent validation or600 biological units.
Raw inputs, gene folds, primary endpoint definitions and family6 bounds remain
fixed; old outputs and the41/100 state failure remain frozen. The panel audit
recounts all public outcomes and distinguishes error suppression from retained
two-clone recovery. All12 inference jobs and audit23121811 pass. The audit checks
600 outcomes,1800 inference-input records,900 unique input hashes and600 exact
gene-fold files. All600 historical fit manifests and10,800 historical output
files match their frozen hashes.

| Condition | Corrected raw units /100 | Corrected effective units /100 | Frozen raw candidate /100 |
|---|---:|---:|---:|
| Matched-null extra profile |0|0|0|
| State-null extra profile |12|12|41|
| Exactly two profiles and exclusively chr12 differences |99|98|98|

The descriptive family6 simultaneous upper bound is4.6747% for either matched
null and21.8909% for either state null. Recovery lower bounds are93.3491% for raw
and91.6455% for effective units. These use the original exact binomial family and
agree with independently inverted binomial tails. They describe exposed
simulations and cannot certify a procedure developed using them. The state
criterion still fails; the unit-invariant arm is not discarded to obtain the
higher raw-arm recovery. Mean all-cell two-clone assignment remains approximately
72.8% (raw0.7280667, effective0.7279833). Detection is not accurate membership.
The output criterion remains exactly the original two-profile/chr12 endpoint;
it does not establish breakpoint or effect-size accuracy.

Authoritative panel evidence is `panel_audit_v2/report.json`. The superseded
pending audit23121770 was cancelled before execution to add frozen-output hash
checks. Two observation-audit failures and two pending cancellations remain in
the task journal, alongside all successful runs; none is hidden or reclassified
as a scientific pass.

No interval-coverage or membership-calibration failure is repaired by this change.
Full-count nuisance fitting, selected memberships and state/model inadequacy
remain open. Biological n=0; independent corresponding DNA and suitable null
replication remain unavailable requirements.
