# Variance decomposition reaches the bulk profile but does not improve it

[Observed pattern, not a guarantee] Composing the existing normal-noise
subtraction with exact posterior boundary averaging gives mean Pearson
0.7648144996 versus the retained additive model's0.7648167082. The fixed primary
fails its point continuation rule. Extended arithmetic confirms the sign and
reproduces every frozen float32 validation window exactly. Reject expansion;
no variance scale, threshold, prior, dataset or alternative estimator is tuned.
All7 CPU jobs complete0:0. Production, old tests and frozen results stay unchanged.

## The implemented intervention and what it tests

The prior normal sampling-variance subtraction changed1131 reference variances
and48 fitted tau values, but hard-path exports were identical. The previous
posterior average instead used additive reference/query variance. Test their
composition to determine whether that export discarded useful variance
information. This is an implemented performance candidate, not a new audit of
unchanged code. [Design frozen before outcomes](../dev/decomposed_boundary_average_20261110/DESIGN.md).

Use all48 saved tumour contexts and20 normal contexts. Reconstruct every stored
decomposition field exactly, retaining original effects, reference gene weights,
beta/tau, state priors, transitions, genomic resets and all interval summaries.
Require total variance=intrinsic reference variance+query sampling variance+tau²
for every query. Average descriptive gene summaries over ALL complete HMM runs,
without truncating posterior mass. Reference-weighted median is the fixed primary;
uniform median and reference mean are ablations, never substitute primaries.

Estimand: regional natural-log RNA abundance relative to the specified normal
reference. Observation model: robust Gaussian bin summaries plus an approximate
reference variance decomposition. Identifying assumptions include adequate
reference composition and working residual variance. Regional transcription,
shared-reference dependence and biased finite-normal variance subtraction remain
alternatives to a genomic explanation. Independent corresponding DNA at the
same resolution distinguishes genomic prediction; technical replication and
repeated-workflow coverage would test the variance interpretation. These
existing patients are regression-only and add no biological replication.

The correction is falsified as an accuracy improvement if it does not beat the
retained primary on fixed paired patients. No mathematical elaboration changes
that rule. The RNA job locks predictions before the separate DNA scorer runs.

## Exact controls, independent numerics and the complete result

The additive replay reproduces all3 previous RNA arrays exactly. Reconstructed
variance components, reference identities and weights are exact. Every one of
the42 earlier method score tables reproduces exactly, retaining their failures.
Keep45 methods, the original24/24 calibration/evaluation patients,23 informative
correlations,13 loss recalls, common genomic support and unchanged scorer.

Independent3-bin/3-state enumeration uses all27 full SciPy Gaussian-density
paths and direct interval medians under both variance choices. Maximum summary
errors are2.78e-17 and4.16e-17 against operation-derived budgets2.81e-11 and
2.74e-11. The variance intervention changes the toy averaged profile by
0.000948701; its zero-normal-noise limit is exact. Existing decomposition,
Gaussian path, median, permutation and reset checks also pass.
[Numerical evidence](../dev/decomposed_boundary_average_20261110/numerics_v1/report.json).

| Method | Mean Pearson /23 | Loss recall /13 | Neutral loss FPR /23 | Precision /calling patients |
|---|---:|---:|---:|---:|
|Retained additive reference median|0.7648167082|0.6027473|0.0437372|0.4837963 /18|
|Primary decomposed reference median|0.7648144996|0.6027473|0.0437372|0.4837963 /18|
|Decomposed reference mean|0.7603644110|0.6137363|0.0537841|0.4467593 /18|
|Decomposed uniform median|0.7513195066|0.6137363|0.0412655|0.3707454 /19|
|Strong pseudobulk, own partition|0.7302613671|0.3996795|0.0412445|0.3585714 /14|

The primary beats41 of42 previous methods in paired mean Pearson, but loses
against the retained best by2.20861124e-6,12 patient wins/11 losses. That fails
the frozen ALL42 requirement. Precision means have different caller populations;
this table does not establish a paired precision benefit. All per-patient values,
including calibration scores, remain in
[patient scores](../dev/decomposed_boundary_average_20261110/score_v1/patient_scores.tsv)
and [full45-method report](../dev/decomposed_boundary_average_20261110/score_v1/report.json).

## Resolve the small difference without changing the model

Original block rounding diagnostics range4.1e-11 to1.8e-6; these are not direct
Pearson bounds. A separately frozen precision follow-up recomputes all48
profiles with63-bit longdouble mantissas, independent direct interval emission
accumulation and cancellation of state-independent Gaussian constants. The
original binary64 model inputs, transition entries and exact dyadic medians
remain fixed. All27 exact Fraction path weights agree within1.08e-19 versus
the derived1.29e-14 allowance.

Maximum profile changes are2.47e-13 additive and3.41e-13 decomposed. EVERY
float32 validation window remains EXACT. Without float32 window rounding,
extended-precision paired mean difference is-2.24116618e-6, the same sign.
Thus the failed point gate is numerically stable. This is not evidence of a
scientifically meaningful tiny decline or validated uncertainty.
[Precision design](../dev/decomposed_boundary_average_20261110/PRECISION_DESIGN.md),
[precision outcome](../dev/decomposed_boundary_average_20261110/precision_compare_v1/report.json).

No reference/patient uncertainty expansion or default adoption follows. The
earlier69-family failure stays unchanged. Clone recovery remains2/5 in
production; these bulk computations provide no clone or spatial validation.

## Every execution and preserved artifact

Research code:
[decomposed_boundary_average.py](../benchmarks/hardening/decomposed_boundary_average.py),
[boundary_precision.py](../benchmarks/hardening/boundary_precision.py).
Source_v1 and source_v2 are immutable; latest validated noneditable wheel,
neutral directories, protected frozen paths, exact code/input/environment/native
hashes and resource records are retained. No dependency or shared install changes.

| Job | Role | Account/QOS | CPU / RAM / requested time | Actual time |
|---|---|---|---|---|
|23157105|independent numerics|nslab|1 /4GB /1h|23s|
|23157106|exact additive replay|innovation|1 /8GB /2h|123s|
|23157107|new composed RNA profiles|nslab|1 /8GB /2h|124s|
|23157108|dependent45-method scoring|nslab|1 /8GB /1h|12s|
|23157230|extended additive replay|innovation|1 /4GB /1h|78s|
|23157231|extended composed replay|nslab|1 /4GB /1h|85s|
|23157232|dependent precision comparison|nslab|1 /4GB /1h|7s|

All CPU/v6&amd, all0:0, no arrays or interactiveQOS. Peak scheduler RSS159436K;
peak actual process RSS218952KB from `/usr/bin/time -v`.
[Exact commands](../dev/decomposed_boundary_average_20261110/jobs.json),
[accounting](../dev/decomposed_boundary_average_20261110/accounting_final.tsv),
[resources](../dev/decomposed_boundary_average_20261110/resource_summary_final.json),
[preservation](../dev/decomposed_boundary_average_20261110/preservation_final.json).
The existing456 Python/4 Rust suite was not rerun for these research-only edits;
the installed wheel identity is verified by every launcher. No new GPU test is
claimed. The full scientific goal remains unmet.
