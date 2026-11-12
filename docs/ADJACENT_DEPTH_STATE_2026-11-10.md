# Composing the fixed corrections fails recovery and specificity

[Observed pattern, not a guarantee] Applying the unchanged depth-conditioned
predictive safeguard to adjacent-reference count folds retains3/5 original
recoveries, rejects2/5, and favors a mixture in one known state-null false
split. Another false split remains numerically unresolved. The frozen screen
fails, so reject composition, threshold tuning and full-workflow expansion.
All9 CPU jobs complete0:0. No production, criterion or frozen-result changes.

## Why this composition was tested

The [depth-conditioned state safeguard](DEPTH_STATE_PREDICTIVE_2026-11-10.md)
with old reference noise rejects all4 prior false splits but fails original109.
The [adjacent-reference estimator](ADJACENT_REFERENCE_NOISE_2026-11-10.md) alone
preserves5/5 original recoveries but fails its fresh phase7 test,3/10 state-null
errors and9/10 power versus production0/10 and10/10. Both failures stay frozen.

Changing the count-size estimator changes the split observations and discovery
q, which could change the guard's power/specificity. This experiment composes
those two independently motivated corrections WITHOUT modifying either: same
estimator, stratum means, discovery-only prevalence, shared contrast/continuous
variance, exact boundary optimizer, starts,24/48/96 quadrature, masks, region
extents and numerical sign diagnostic. It is not a chosen random redraw or a
noise value selected against truth.

Before any composed scores, lock ALL5 new original recovery cores and ALL3
new phase7 false-split cases, exact identifiers in CASES.json. Each has exactly
one accepted chr12 event/common-support confirmation core. Model fitting uses
discovery counts only; parameters lock before confirmation entries decode.
Full=D+C holds exactly. [Frozen design](../dev/adjacent_depth_state_20261110/DESIGN.md),
[all8 case paths](../dev/adjacent_depth_state_20261110/CASES.json).

The estimand and identifying limitations are unchanged: conditional confirmation
RNA predictive advantage, with learned sizes/prevalence and fixed upstream
selection. Continuous transcription/reference variation can be RNA-equivalent
to a genotype contrast. This selected exposed panel has biological n=0 and
cannot provide a5% error estimate, DNA validation or cell membership truth.

## Full fixed-panel outcome

| Case | Cells | Genes | Mixture predictive advantage | Numerical diagnostic | Safeguard |
|---|---:|---:|---:|---:|---|
|original_20261031|508|100|4.381734208|5.10984838e-05|retain|
|original_20261107|512|93|4.856707025|7.06807e-05|retain|
|original_20261108|514|94|6.772815641|7.17585137e-05|retain|
|original_20261109|502|93|-0.209526449|4.05057327e-05|reject|
|original_20261110|510|94|-0.708815903|5.36372306e-05|reject|
|state_202611097002|600|93|-2.383041649|0.000129602969|reject|
|state_202611097005|600|93|1.158607826|9.27987064e-05|retain|
|state_202611097007|598|94|-0.000000528|9.47717456e-05|unresolved|

Every original must have resolved-positive advantage and every false split
resolved-negative before expansion.109 is-0.209526449 and110 is-0.708815903,
while false split202611097005 is+1.158607826. No ordinary increasing cutoff on
THIS score can retain all originals and reject that false split. This is an
ordering observation on the fixed data, not a proof that another method or
independent genomic measurement cannot distinguish them.

State202611097007 is-5.28e-7 versus diagnostic9.48e-5; explicitly unresolved.
The diagnostic is a successive-order convergence heuristic plus rounding, not
a formal quadrature enclosure or biological interval. No higher precision run
is needed for the expansion decision: the other three decisive wrong-sign
cases already fail by margins far beyond their diagnostics.
[All fits, order-specific predictions and decisions](../dev/adjacent_depth_state_20261110/summary_v1/report.json).

No full guarded workflow or fresh phase8 calibration is run after this failure.
The rejected unguarded phase7 final100 continuation is not reopened. Original
recovery assertions and one-sided simultaneous95% whole-specimen error upper
<=5%, with power and abstention, remain unchanged. Production recovery is2/5;
no tested research method yet satisfies the joint recovery/error objective.

## Numerical evidence, execution and preservation

The model source is BYTE-IDENTICAL to the independently validated
[depth_state_predictive.py](../benchmarks/hardening/depth_state_predictive.py).
Its exact Fraction probabilities, independent full-PMF adaptive quadrature,
heat-equation gradients, nested null and permutations are reused explicitly
through a read-only numerical-evidence symlink. Those checks are NOT rerun
merely for a new data panel. Input manifests resolve/hash the existing evidence.
There are no new model constants, production dependencies, old test edits or
shared installations. The private wheel remains the one validated by456 Python/
4 Rust tests; no new package test suite is claimed this turn.

All CPU/v6&amd:8 workers1CPU/8GB/2h alternating innovation/nslab and1 dependent
summary nslab1CPU/4GB/1h, no arrays/interactive QOS. Immutable source_v1, neutral
working directories, protected frozen paths, wheel/native/source/input hashes,
actual backend and resources remain.

| Job | Role | Account/QOS | Runtime | Outcome |
|---|---|---|---|---|
|23156948|case_00|innovation|00:01:33|complete0:0|
|23156949|case_01|nslab|00:01:22|complete0:0|
|23156950|case_02|innovation|00:01:09|complete0:0|
|23156951|case_03|nslab|00:01:09|complete0:0|
|23156952|case_04|innovation|00:01:52|complete0:0|
|23156953|case_05|nslab|00:01:54|complete0:0|
|23156954|case_06|innovation|00:01:46|complete0:0|
|23156955|case_07|nslab|00:03:10|complete0:0|
|23156956|summary|nslab|00:00:10|complete0:0|

Peak scheduler RSS 118848K (23156953.batch).
[Exact commands](../dev/adjacent_depth_state_20261110/jobs.json),
[all accounting](../dev/adjacent_depth_state_20261110/accounting_final.tsv),
[resource summary](../dev/adjacent_depth_state_20261110/resource_summary_final.json).
All original package/native/test files and frozen outputs are preserved;
this is an executed negative method-development result, not a promoted fix.
