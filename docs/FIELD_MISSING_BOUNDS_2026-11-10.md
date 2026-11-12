# Fixed-cohort mean advantages survive the missing comparator scores

Job23167023 completed0:0 in16s. Every deterministic check passes. On the fixed
181-line field table, the augmented RNA stack has a positive mean-score gap
against every one of the13 other recorded paths on both endpoints, including
the worst possible DEFINED scores for missing CopyKAT/SCEVAN populations.
This is a sensitivity result on an exposed dataset, not a new significance
test, independent superiority finding, or restoration of missing predictions.

| Comparator | Known lines | Lowest possible mean1Mb Pearson gap | Lowest possible mean arm Spearman gap |
|---|---:|---:|---:|
|CopyKAT|168|+.0705627|+.0978198|
|SCEVAN|167|+.1063745|+.1451828|
|inferCNV R|181|+.0393384|+.0373573|
|infercnvpy|181|+.0834843|+.0520468|
|Frozen HaploDose gene-ratio path|181|+.1269786|+.0124255|
|Smoothed pseudobulk5|181|+.0445580|+.0143231|
|Gaussian-MDL pseudobulk|181|+.0404409|+.0221173|
|Gene-level pseudobulk|181|+.3552363|+.0021384|
|Previous RNA stack|181|+.0041340|+.0024403|

The complete report retains all26 contrasts, all181 cases, all14 methods and
every missing identifier. The primary augmented profile has a finite score
for every case. For missing comparator scores, the mean-gap lower endpoint
assigns correlation+1 and the upper assigns-1; all other scores remain exactly
as recorded. No profile, label or score is filled in the source table.

All input scores become exact dyadic fractions. Sums and division by181 are
exact, and printed endpoints round outward. Exhaustive endpoint/zero
completions on a small independent fixture reproduce the extrema; complete
identity and directed printing checks pass. The ranges are fixed-dataset
identification ranges, not confidence intervals. The small arm advantage
against gene-level pseudobulk remains small and does not imply a population
effect. Shared references, repeated exposure and uncertain DepMap truth remain.

The original26-test field family STILL fails its complete-output requirement.
No p-value, Holm family, coverage rule or formal verdict changes. These bounds
only show that unknown defined comparator scores cannot reverse these observed
mean rankings. If a method's output remains undefined, it remains unavailable;
the hypothetical finite extension is not a deployable prediction. Future
cohorts and biological validation are still required for a SOTA claim.

No inference, refitting, method selection, default, scientific criterion or
frozen result changes. Biological replication added=0. The job uses CPU/
innovation,1CPU4GB30min,read-only validated environment and911-file snapshot
retaining all909 prior source files exactly.

[Post hoc design](../dev/field_missing_bounds_20261110/DESIGN.md),
[exact calculation](../benchmarks/hardening/field_missing_bounds.py),
[all26 ranges and missing IDs](../dev/field_missing_bounds_20261110/bounds_v1/report.json),
[complete manifest](../dev/field_missing_bounds_20261110/bounds_v1/run_manifest.json),
[submission](../dev/field_missing_bounds_20261110/jobs.json),
[unchanged original comparison](GENE_LAPLACE_PROFILE_2026-11-10.md).
