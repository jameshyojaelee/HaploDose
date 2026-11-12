# Conditional likelihood does not identify common tumour dosage

Correction to the profile-trace interpretation: the historical profile counts
below score raw returned profiles, including the normal row that refinement
actually fixes to zero. They are not counts of applied profile transitions.
Unchanged-production replay finds 31/67 applied decreases versus 22/67 raw
decreases, with all original recovery criteria unchanged. Applied counts for
the other research arms were not run. Historical artifacts remain unchanged.
See [exact transition and independent replay](PROFILE_TRANSITION_2026-11-09.md).

The offset-only optimizer experiment does not rescue clone recovery. It separates
three questions: does the code evaluate its stated objective correctly, do
updates increase that conditional objective, and do they recover the supported
clone structure? None implies the others. No package default is changed.

## Contract and competing explanations

Conditional estimand: shared regional RNA log-rate offset, with responsibilities,
reference, depth, programmes, clone effects and dispersions held fixed. The
working objective is Q=sum_ik R_ik sum_j log NB(Y_ij;m_ikj,phi_j), where
m_ik=L_i softmax(log base+offset+dose_k+depth_i+programme_i). Keep the original
retained weight mass. Fix the arbitrary constant-logit gauge during optimization;
median centring does not identify the biological meaning of a regional effect.

The strongest competing explanation for a common tumour RNA effect is a shared
nuisance offset. With zero reference-component weight, adding a regional vector
to the offset and subtracting it from each tumour profile leaves every tumour
mean identical. No RNA likelihood optimizer distinguishes those labels. Positive
reference-component weights inform this direction only through reference-weighted
observations; fitted reference-like membership is not independent normal truth.
This is a conditional identification counterexample, not evidence that every
real specimen has zero normal support. Independent reference identity and
transport evidence, followed by corresponding DNA, would distinguish biological
dosage from a technical/reference explanation.

## Design and independent calculations

Task `dev/offset_objective_20261109` preserves the current dirty baseline and uses
the preceding validated wheel read-only. No build, dependency, shared-installation
or preexisting test change. Source snapshots and all outputs have new identities.

The research L-BFGS-B solver has an analytic composition-normalized NB gradient,
one fixed gauge coordinate and explicit convergence/objective diagnostics.
Its objective uses clone-wise cell-by-bin arrays; it does not allocate the full
cell-by-clone-by-bin tensor. The objective is exact for its stated conditional
working model; there is no theorem of global convergence in general dimensions,
nor an exact joint-EM claim. Numerical convergence is distinct from recovery.

Twenty-four constructor-compatible integer-count cases compare its two-bin
solutions against independent SciPy NB logpmfs and Brent searches, refined
finite differences, and scalar Decimal80 NB likelihood ratios. Existing physical
exposures, nuisance contrasts, hard/soft/dropped weights and fixed count
perturbations remain. Iteration caps bound resources; they do not certify
stationarity. Tolerances derive from arithmetic bounds and independent numerical
refinement, not a scientific decision threshold.

All five original recovery inputs are run under current production, the rejected
conditional-ratio update and conditional NB optimization. Only the offset differs.
Instrumentation also records each profile update's Q change with responsibilities,
nuisance and dispersion fixed. Save the most negative profile and offset context
for every fitted model, selected without truth. Independently replay complete
contexts using full NB logpmfs and bounded cell subsets using Decimal80. Recount
all six original criteria for all15 arm/case pairs. Diagnostic-only arms must
reproduce the earlier complete criteria, not merely their pass counts.

## Failures retained and numerical correction

Three initial trace jobs completed their inference panels but failed while
serializing a NumPy boolean. Their outputs and failed job records remain. The
recorder correction changes representation only. An audit encountered the same
serialization issue; its corrected replay passes the original90 recovery criteria
and90 worst-context checks. No scientific criterion was relaxed.

The subsequent identification probe exposed cancellation in the first research
Q implementation: subtracting separately normalized log probabilities and
logaddexp terms could miss the original bound when the true likelihood change
was zero. The initial96-case probe fails and is preserved. An explicit archived-
helper replay records that failure; its successful execution is not a passing
numerical result. The corrected helper forms normalized differences and NB
ratios using log1p/expm1. The same96-case probe and original24 optimization checks
then pass, without enlarging their bounds. All three workflow arms are rerun
with the corrected arithmetic. This correction affects research code only.

The identification probe independently verifies768 tumour probability vectors
(1536 components),96 conditional-Q comparisons and linear reference-mass scaling.
At zero anchor mass the common-dose/offset direction is flat; with positive
anchor mass and counts equal to the reference mean, moving from that mean lowers
Q. This proves a limitation of the observation model, not population performance.

## Workflow disposition

Final stable-arithmetic jobs23095491/492/493 and audit23095494 complete. The audit
reconstructs all90 original recovery criteria, replays90 worst contexts using
full NB logpmfs and bounded Decimal80 calculations, and verifies that both
diagnostic-only arms reproduce the earlier complete criteria exactly.

| Arm | Recovery | Profile Q decreases | Offset Q decreases |
|---|---:|---:|---:|
| Current production |2/5|22/67|33/67|
| Conditional count ratio |1/5|26/69|0/69|
| Conditional NB optimizer |1/5|26/63|0/63|

These are update counts from dependent development trajectories, not independent
replicates or error probabilities. The minimum profile Q changes are−35.4967,
−13.1949 and−13.1489 weighted log-likelihood units respectively. NB offset gains
are positive throughout(minimum0.199812), with no reported optimizer failure.
Both conditional arms pass only20261108; current production passes20261031 and
20261107. Thus an offset-step likelihood decrease does not explain the observed
conditional-ratio regression, and improving that step is insufficient to rescue
recovery. This does not exclude a contribution from approximate offsets in other
cases, or prove which remaining stage causes all three original failures.

Stable numerical gate23095489 passes all24 original cases. Identification
gate23095488 passes all96 cases under unchanged bounds. Archived-helper
replay23095487 records all_checks_passed=false as expected. Its zero exit code
means the old failure was reproduced; it is not a passing accuracy result.

Segmentation, shared tumour profiles and nuisance estimation need not maximize
this unpenalized cell likelihood. A decrease is not automatically a coding bug:
segmentation restrictions or priors could justify it. The package would need a
stated joint objective or an independently demonstrated scientific advantage to
justify those tradeoffs. Replacing only the offset with a better conditional
optimizer cannot supply that argument.

Retain the original three recovery failures, frozen false-split/power results,
5% specimen criterion and simple baselines. This experiment establishes neither
calibrated uncertainty nor a biological clone-error guarantee. Biological n=0;
former evaluation datasets remain regression-only. A common tumour RNA effect
cannot be promoted to DNA dosage by choosing a different decomposition between
clone effects and a fitted shared offset.

## Reproducibility and next decision

All17 CPU jobs are terminal:12 complete and5 preserved failures(four recorder
serialization failures, one original numerical counterexample). The final
seven-job chain has no execution failures. `jobs.json` and `accounting_final.tsv`
record resource headers, accounts, runtime, step memory and dependency status.
The final NB process uses1527684KiB peak RSS and82.91s; final audit229592KiB and
36.06s. Both innovation and nslab CPU allocations were used. No GPU or other
partition was needed. A static comparison in `source_integrity.json` verifies
all76 initial package/native/test/README/build files unchanged.

All jobs use the read-only wheel/environment from
`dev/offset_expectation_20261109/validation_23092399`, wheel SHA256
`711325c2020268d339822afc896efbd36c1f1370751317e32368e7eabb55b93f`.
Exact source, native, wheel, environment and argv records are in each
`analysis_JOB/execution.json`. Final source is `stable_source_v1`; numerical,
identification, workflow and audit outputs have distinct directory identities.

The rejected optimizer is retained as research code. Full suite rerunning is
not needed for an unchanged package; the actual frozen-input workflows and
independent numerical tests are the checks for this research addition. No new
null/power performance claim is made for a candidate that already fails the
unchanged primary recovery contract. Earlier600-case false-split/power results
and168-bound family remain unchanged.

The next decisive comparison is a coherent reference-anchored model against
the simpler restriction of no free regional shared offset, with explicit
reference-transport assumptions. Test common-dose preservation and clone
contrasts together. Any improvement must survive the original recovery and
complete-workflow null/power gates; better conditional Q is already shown to
be an inadequate promotion rule. Independently justified normal/reference
measurements, rather than fitting a convenient offset, must supply the anchor.
