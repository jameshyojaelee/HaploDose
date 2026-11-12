# Removing approximate normal sampling noise does not improve the profiles

The implemented variance decomposition changes all1131 normal-bin variances and
all48 fitted query extra-variance values, but all three exported RNA profile
arrays equal the preceding Gaussian model exactly. Mean weighted-median Pearson
remains.760619. The declared positive-gain prerequisite fails; no larger
reference/patient uncertainty experiment follows. Retain this negative result
without changing production, tests, thresholds or frozen artifacts.

## Scientific model and decisive intervention

The original Gaussian boundary model uses empirical normal-bin spread PLUS
query sampling variance. Its design explicitly warns that normal sampling noise
may be counted twice. The
[new design](../dev/normal_noise_decomposition_20261109/DESIGN.md) tests that
specific mechanism while retaining all other parameters and the33 preceding
comparators. This is a scientific approximation, not a computational bug fix.

Estimand remains natural-log regional RNA abundance relative to the declared
normal reference. The working observations are robust gene summaries with a
Gaussian bin model. For each normal, recompute its existing leave-one-out robust
bin summary and approximate weighted-quantile sampling variance. Use the
bin-wise median across20 normals as typical normal sampling variance. Define
intrinsic spread as max(empirical normal-bin MAD variance - typical sampling
variance,0), then add the query sampling variance. The zero boundary expresses
the variance parameter space; it is not an epsilon. Fit the SAME adjacent-bin
extra variance and SAME beta grid, priors, transitions and genomic resets.

The identifying assumptions are approximately Gaussian log-expression,
representative normal measurement noise and an adequate decomposition despite
heteroscedasticity, finite-n MAD bias, shared references and correlated genes.
Subtracting approximate measurement variance from robust spread is not an
unbiased variance-component estimator. Reference transport and transcriptional
programmes remain alternatives to DNA dosage. Corresponding independent DNA at
common support tests profile prediction; this exposed panel cannot establish
independent biological validation or interval coverage. Biological n=0 is added.

[Implementation](../benchmarks/hardening/normal_noise_decomposition.py) is
research-only. Keep the original zero intercept, genes, reference roles,
corrected effects, weights,48 tumours,20 normal rows,24/24 calibration/comparison
patients and failure/undefined rules. No DNA outcome chooses a reference or
parameter. All RNA predictions and model decisions freeze before DNA scoring.

## Computation and actual performance

Before new primary fitting, all48 original fixed-zero Gaussian contexts and
their beta/tau/grid/near-tie/rounding/calibration decisions reproduce EXACTLY.
All33 preceding method scores reproduce exactly;36 methods total. Independent
binary-rational variance identities, no-noise limit, nonnegative boundary,
normal-summary identity, gene permutation and original complete-path Gaussian
likelihood checks pass. Numerical checks do not claim sampling calibration.

| New output | Mean Pearson,23 patients | Gain versus strongest pseudobulk | Gain versus preceding Gaussian weighted median |
|---|---:|---:|---:|
|Uniform median|.744107|+.013846|-.016512|
|Reference-weighted median|.760619|+.030358|0|
|Reference-weighted mean|.755013|+.024752|-.005606|

Weighted-median loss recall stays.540476 on13 defined patients; FPR stays.042373
on23. Its gain versus the strongest pseudobulk remains a point estimate, not a
new superiority result. The earlier42-family combined interval/Holm failure
remains. No candidate gains over ALL native, strongest pseudobulk, normal-centered,
original Gaussian and median-intercept Gaussian baselines, so reject continuation.

[Primary comparisons and checks](../dev/normal_noise_decomposition_20261109/bulk_v1/report.json).

Identical scores alone could indicate a wiring error. A separately declared
[stored-output check](../dev/normal_noise_decomposition_20261109/EFFECT_DESIGN.md)
therefore verifies the intervention without rerunning inference. All1131 normal
intrinsic variance estimates change, none clips to0, and maximum removed
variance is.0007923961720693352. All48 query tau values change; no fitted beta
changes. All three complete exported RNA arrays, including normal-row NaNs,
equal their preceding outputs EXACTLY. The intervention reached the model but
does not improve these predictions. Do not substitute changed internal
posteriors for an accuracy gain.
[All stored variance/path/posterior/export comparisons](../dev/normal_noise_decomposition_20261109/effect_v2/report.json).

## All execution and preservation

| Job | Role | Account/QOS; CPUs / RAM / requested time | Outcome |
|---|---|---|---|
|23152892|independent numerical checks|nslab;1 /4GB /1h|complete0:0,15s|
|23152893|exact controls and primary panel|innovation;1 /8GB /1h|complete0:0,100s; gain gate fails|
|23153071|stored intervention verification|nslab;1 /4GB /1h|complete0:0,6s; exports exact|

Source_v1 runs the first two jobs. Source_v2 adds only the stored-output
verification. The existing validated noneditable wheel/native binary is read
only from neutral compute working directories. Exact commands, allocations,
backend, versions, process resources and immutable sources live in the
[task directory](../dev/normal_noise_decomposition_20261109/).
No production/native/original-test/dependency/shared-installation/frozen writes.
Full package/GPU tests, broader uncertainty and biological validation are NOT
RUN. The separate fixed-platform clone experiment remains queued on its
historical node; this profile result does not repair clone recovery or the5%
whole-specimen false-split/power failure.
