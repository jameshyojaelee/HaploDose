# Matching the robust intercept does not improve the retained profile

The median-intercept Gaussian boundary candidate reaches mean Pearson.757306,
below the preceding fixed-zero model.760619. All48 previous Gaussian contexts,
their beta/tau/grid decisions and all30 preceding method scores reproduce
exactly. Reject the new candidate at its unchanged point-gain prerequisite;
larger reference resampling is not run. Higher loss recall accompanies higher
FPR. No default,production,original test,dependency or frozen change.

## Scientific mechanism and fixed correction

Gene effects are centered by UNWEIGHTED within-bin medians, while the boundary
HMM observes WEIGHTED bin medians. Those locations need not agree. The
[locked design](../dev/centered_gaussian_regions_20261109/DESIGN.md) changes
only the boundary intercept from0 to the median of the actual robust bin
summaries. It adds no smoothing/effect grid and does not recenter exported
gene amplitudes. This is a robust location estimate, NOT the joint Gaussian
mixture MLE. The old fixed-zero path is retained exactly.

Estimand remains regional natural-log relative RNA abundance against the
declared normal gene reference. The working Gaussian observation model and
reference/count variance approximations are unchanged. Majority reference-like
RNA regions identifies this median intercept. Genome-wide transcription,
reference composition and purity/ploidy ambiguity are alternatives. Absolute
genome doubling remains unidentifiable. Corresponding independent DNA at common
support tests prediction; biological coverage and replication require additional
evidence. This is exposed development/regression evidence, not new evaluation.

Keep all prior normal gene effects,weights,LOO variability,query variance/tau,
state levels,priors,p_change,beta grid,genes,genomic resets and sample roles.
Three fixed outputs use the same new regions: uniform median,weighted median,
weighted mean. The same48 tumour patients,20 normal rows and24 calibration/
24 comparison split remain, with23 informative correlations/13 loss recalls.

## Independent computation and exact preceding controls

[Intercept implementation](../benchmarks/hardening/centered_gaussian_regions.py)
uses an explicit optional mode in the research
[Gaussian fitter](../benchmarks/hardening/reference_gaussian_regions.py).
Independent scalar median checks and an exact constant-shift fixture verify
the intercept,posteriors,MAP path and all beta scores. The preceding fixed-zero
model fails that invariance fixture. This is a computational property, not proof
that the new intercept better predicts DNA. Existing independent complete-path,
variance,normal-LOO,Fraction and gap/permutation checks pass.

Before new primary fits, recompute ALL48 preceding fixed-zero Gaussian contexts
and every stored beta,tau,grid score,near-tie,rounding-budget field EXACTLY.
Then save all new RNA predictions/decisions before loading DNA scores. All30
preceding method scores are exact;33 methods total. Keep all patients,callable
support,failures and undefined endpoints. No missing comparator can improve a
verdict. [Exact-control record](../dev/centered_gaussian_regions_20261109/bulk_v1/fixed_zero_control.json).

| Method | Pearson(23) | Loss recall(13) | Classified precision | Loss FPR(23) |
|---|---:|---:|---:|---:|
| Median intercept,uniform median |.741028|.573443|.332744(18)|.044675|
| Median intercept,weighted median |.757306|.650366|.402024(19)|.050221|
| Median intercept,weighted mean |.755572|.624725|.425000(18)|.048541|
| Fixed-zero Gaussian weighted median |.760619|.540476|.384804(17)|.042373|
| Strong pseudobulk / own partition |.730261|.399679|.358571(14)|.041244|

Median-intercept weighted median has10 correlation wins/8 losses/5 ties versus
the preceding Gaussian model, but the paired mean is-.003313. Its gain versus
strong pseudobulk is+.027045(18 wins/5 losses), insufficient for the declared
continuation rule requiring a positive mean against ALL preceding Gaussian,
normal-centered,native and strong pseudobulk baselines. No candidate passes.
Do not change the failed correlation prerequisite after seeing higher recall.
Precision is conditional on making calls and has different patient populations;
it cannot establish unconditional accuracy. Achieved FPR is descriptive, not a
5% specimen or biological guarantee.

[All33 methods and primary contrasts](../dev/centered_gaussian_regions_20261109/bulk_v1/report.json).
Reject as a default replacement; retain the lower-correlation/higher-recall
tradeoff as a negative result. Reference resampling,bootstrap superiority,
interval coverage and new biological validation are explicitly not run after
the failed prerequisite. The earlier Gaussian42-family failure remains.

## Compute and preservation

| Job | Role | Account/QOS; CPUs / RAM / requested time | Outcome |
|---|---|---|---|
|23151881|independent numerical checks|nslab;1 /4GB /1h|complete0:0,20s|
|23151882|exact controls and primary panel|innovation;1 /8GB /1h|complete0:0,100s; point gate fails|

CPU compute node006,neutral cwd,existing read-only noneditable validated wheel
and native binary. Peak process RSS249216KiB. Shared pre-edit baseline and
source_v1 live in the predictive-gene task; the
[profile task](../dev/centered_gaussian_regions_20261109/) records exact command
arrays,backend,environment,resources and hashes. No arrays,GPU,installation,
native/production/original-test edits or frozen writes. Full package tests are
not rerun for research-only code. Clone recovery,biological validation and the
5% whole-specimen false-split/power objective remain unresolved.
