# Adjacent-depth noise experiment stops at its exact recovery control

The new adjacent-depth reference noise estimator passes its independent numerical
identities, but its full historical recovery control fails. Both dependent
interventions are cancelled before execution. Their recovery and false-split/power
behaviour are not measured. Do not describe them as scientifically rejected or
as fixes. The five control runs retain the same 30 criterion records and 4/5
research recovery, but small output differences prevent the stronger exact
reproducibility claim. Production remains 2/5; all three original production
failures and the 5% specimen error/power objective remain unresolved.

## Intended performance intervention and falsifier

The existing reference-quintile estimator gives common gene NB sizes 8.36–8.99.
Residual depth-related composition within a quintile could inflate its noise
estimate and remove membership information after count splitting. This is a
specific competing explanation, not an established cause of failed recovery.
The [frozen design](../dev/paired_reference_noise_20261109/DESIGN.md) replaces only
that estimator with disjoint adjacent-depth pairs within dataset/reference type,
breaking ties by literal cell ID. One terminal triplet retains every cell when
the population is odd. No arbitrary bandwidth, query count, holdout reference,
DNA label or true dispersion chooses pairs or sizes.

For fixed exposures L and reference counts Y in a pair/triplet, compute
A=sum Y(Y-1) and B=sum L^2 times the cross-cell gene-product sum divided by the
cross-cell exposure-product sum. Pool A-B and B across blocks and set
phi=sum(A-B)/sum(B), size=1/phi. Long-double accumulation reduces cancellation.
Unsupported cross-cell mean support or nonpositive dispersion fails explicitly;
there is no clipping or oracle fallback. The implemented routines are
[paired_reference_noise.py](../benchmarks/hardening/paired_reference_noise.py)
and [paired_reference_recovery.py](../benchmarks/hardening/paired_reference_recovery.py).

The estimand is common reference-gene NB dispersion transported to split-count
likelihoods. Under independent gene NB observations, fixed exposure and constant
composition within each block, E[A]=(1+phi)E[B]. Adjacent observed depth does not
prove constant composition, and retained-gene totals are random rather than
exogenous exposures. Reference programmes and transport failure remain competing
explanations. Independent RNA replication can examine that approximation;
corresponding DNA and biological replication are needed for genomic clone claims.
No biological specimen is added here.

Independent Fraction cross-cell moments, the conditional NB expectation identity,
reference-order invariance, odd-triplet retention and unsupported-reference
rejection pass. The operation-count/magnitude rounding allowance is
3.2250500458e-14, not a new scientific decision tolerance. These checks establish
the stated formula, not its unbiasedness with observed totals or its power.

Had the control passed, the two planned interventions would have used this noise
estimate with unchanged weighted confirmation/merging, and with the already fixed
count-gene shared-membership confirmation. Both require 5/5 original recovery
before a larger fresh null/power experiment. The two-candidate/two-null/two-look
family would use alpha .00625 per exact one-sided bound, requiring an upper
specimen false-split bound at most 5%, alongside power/membership/abstention.
No posterior support, matching, merge veto, initial-only acceptance, reference
selection, original assertion or stopping rule is relaxed.

## Scheduling amendments precede the first and only control execution

The initial request is 4 CPU/16GB/3h on the previously successful ne1dc6-006 node.
That node is fully occupied with a next-morning start estimate. While all three
jobs remain unstarted, shorten time to 30min, still over six times preceding
runtime. The recorded pending allocation first moves to idle 001, which is then
found reserved for gglab_cpu and unavailable. Before any execution, move to 010
with 2 CPU/16GB/30min. Previous controls used about 99% and 121% CPU. No algorithm,
seed, solver tolerance or exact output check changes. Allocation controls follow
the official [scontrol interface](https://slurm.schedmd.com/scontrol.html).

The [allocation amendment](../dev/paired_reference_noise_20261109/NODE_ALLOCATION_AMENDMENT.md)
explicitly limits this experiment to one actual control attempt: failure stops
the candidate chain without a search across nodes for a passing result. All
pending-only resource/node/CPU update command records remain saved. The amendment
changes the original node/thread allocation for availability, so it is not a
literal original-allocation replay. Actual inference executes only on 010.

| Job | Role | Final resources | Outcome |
|---|---|---|---|
|23149990|independent numerical gate|nslab CPU1/4GB/1h|complete 0:0,20s|
|23149991|full historical control|innovation CPU2/16GB/30min,010|failed 1:0,99s|
|23149992|new size, weighted stages|nslab CPU2/16GB/30min,010|cancelled, never executed|
|23149993|new size, count-gene stages|innovation CPU2/16GB/30min,010|cancelled, never executed|

## What actually agrees and what fails

All five native inferences finish. All 30 original scientific criterion records
agree exactly with the historical reference-quintile research control. Recovery
is still 4/5 and 109 still fails. The exact table comparison stops on the first
seed's log.tsv: historical 196 rows versus new 197 rows. The failed exception,
nested results and dependent-job cancellation records remain intact.

Bounded file/text inspection after failure finds 72/75 byte-identical scientific
arrays/tables. Three tables in seed108 differ: arms.tsv has three changed numeric
fields, profiles.tsv.gz has 43, and segments.tsv has 11. Example mean RNA dose
.08193276077508926 becomes .08193273097276688. These small differences do not
change the 30 recorded scientific criteria, but they invalidate exact output
identity. Some log differences remain after removing numerical-ambiguity actions
in seeds108/110. Do not import the previous experiment's stronger "only logs
differ" characterization into this run. Cause is not isolated, and neither
threads, hardware nor the coincident-mixture identity below establishes it.

[Failure characterization](../dev/paired_reference_noise_20261109/control_failure_characterization.json)
records every file and bounded example. Removing ambiguous rows is diagnostic
inspection only; it does not replace or pass the frozen exact-log gate. No new
candidate dispersion, recovery estimate or false-split/power rate exists. Both
interventions remain **not run**, rather than assigned an invented failure rate.

## Separate bounded mixture probe, without a production correction

The [probe design](../dev/mixture_coincident_20261109/DESIGN.md) tests a smaller
identity: a mixture of two identical null components has null-relative objective
zero, mixing-logit derivative zero and posterior membership equal to its prior.
Use 17 cells, two features, O=E=2, NB size2 and seven predeclared logits. Fraction
weights independently sum to one. A log1p/expm1 and logistic-difference prototype
is evaluated alongside the existing implementation.

Job23150223 (nslab CPU1/4GB/30min) completes 0:0 in16s. All seven existing objective
values are exactly zero. Four configurations show floating membership/gradient
identity differences: maximum gradient1.8873791419e-15 and prior discrepancy
1.1102230246e-16. The prototype returns exact zero/identical priors for this
fixture. The existing operation-derived gain budget already rejects tiny
unresolved scientific gains. This result does not demonstrate lost recovery,
isolate the historical control discrepancy or justify changing production merely
to remove a diagnostic row. No core numerical change is made.

[Probe script](../benchmarks/hardening/mixture_coincident_probe.py),
[full results](../dev/mixture_coincident_20261109/probe_v1/report.json).
A future numerical correction would need a consequential counterexample,
independent likelihood/gradient agreement, the unchanged tests, an isolated new
wheel and original recovery/error-power checks. The present failed gate remains
failed. No node retry, baseline rewrite or relaxed tolerance is performed.

## Artifact identity, scientific status and next constraint

Both tasks use the validated read-only native-region wheel and neutral-cwd import
checks. The paired task executes immutable source_v1. The probe executes the
count-weighting source_v1 snapshot, explicitly recorded in jobs.json, because it
contains the same unchanged installed Python/native sources and new probe script.
All designs, scripts, command arrays, accounting, process memory/runtime, backend,
snapshot integrity and initial dirty state remain in the task folders. The
earlier 421 Python/4 Rust passes and one GPU-only skip are separate evidence;
full suites are not rerun for these research-only modules. No shared installation,
production/native/original test, frozen output, scientific criterion or dependency
changes.

Computational status: moment identities pass; full exact control fails; candidate
execution stops. Identification status: conditional working noise transport is
unestablished, and RNA-defined populations do not identify genomic clones against
equivalent transcriptional states. Biological status: no new corresponding DNA,
independent null specimens or unopened validation panel. The earlier 4/5 candidate
remains rejected for 13/100 state-null false splits, and production remains 2/5.
An unchanged scientific decision rule, suitable independent truth and simultaneous
error/power assessment must establish a fix; more CPU alone cannot do so.
