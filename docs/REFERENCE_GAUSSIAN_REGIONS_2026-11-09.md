# Robust RNA boundaries improve a point profile but do not pass superiority

The new reference-scaled Gaussian boundary model reaches mean comparison
Pearson.760619, versus.756351 preceding centered weighted median,.738056
installed native median and.730261 strongest pseudobulk. Its42-contrast
uncertainty analysis does not pass the unchanged combined interval/Holm rule
against any of those three. Retain experimental results; do not promote it.
Loss recall decreases relative to the preceding candidate. All12 CPU jobs
complete0:0 and the full30-method original-reference replay is exact.

## Scientific question and implemented model

The preceding estimator uses robust gene effects but boundaries learned from
summed-gene counts. A few highly expressed genes can dominate those boundaries.
The [locked primary design](../dev/reference_gaussian_regions_20261109/DESIGN.md)
changes boundary learning while keeping gene effects, normal centering,
reference/count weights, aggregation support and genomic resets fixed.

Estimand: regional natural-log RNA abundance relative to typical normals.
Observation model: continuous bin contrasts with a Gaussian working error
around the existing four beta-scaled RNA-response levels. This does not treat
gene medians as raw molecule counts. Normal reference transport, sparse
breakpoints and locally stationary, independent bin errors are assumptions.
Correlated regional transcription, tissue composition and reference mismatch
can instead explain the observations. Corresponding independent DNA at common
resolution tests prediction, not absolute DNA copy number or cell clone identity.

For each bin, compute the existing reliability-weighted median of its centered
gene effects. Estimate bin spread from20 conditional normal profiles: each
target normal's gene location AND weights use the other19 normal effects.
Those inputs retain their original leave-one-normal nuisance fits and share
other reference data; they are not independent confirmations. Spread is
squared MAD/Phi^-1(.75), without subtracting normal measurement noise.

Add the independent-Gaussian-log-count asymptotic weighted-quantile variance
pi/2*sum(a_g^2)/(sum(a_g*sqrt(E_g+.5)))^2, with the unchanged gene weights
a_g=1/sqrt(v_g). This is a working approximation derived from weighted-quantile
densities at a common location, not exact count uncertainty. Small gene bins,
overdispersion, correlation and possible double-counted normal noise limit it.
Additional query mismatch comes from adjacent-bin differences after subtracting
the median adjacent working variance, with the existing tau_floor.02.

The defined gene contrast already fixes centre0. Select beta by Gaussian HMM
marginal likelihood on the SAME existing bulk grid/bounds; retain native state
levels, normalized prior,p_change and chromosome/arm/gap resets. Float64
forward/backward and Viterbi record every grid score, selected beta/tau and
numerical near-ties. State weights remain unvalidated working support. No
calibrated confidence, DNA-state or detection-limit claim follows.

[Implementation](../benchmarks/hardening/reference_gaussian_regions.py)
feeds the research reference helper behind an explicit Gaussian-region option.
Export uniform median, reference-weighted median and reference-weighted mean
on the SAME new regions. Uniform aggregation is a mandatory matched baseline.
Default native/reference/count modes remain unchanged; no production/native
or original-test edits, new dependency or frozen-output changes.

## Independent checks and complete point results

Independent SciPy Normal densities plus complete4^3 path enumeration verify
evidence, posteriors and MAP. Maximum posterior discrepancy2.22e-16 fits an
operation-derived allowance8.12e-12. Separate scalar quantile-variance and
LOO-normal fixtures, gene permutation and exact genomic-block decomposition
pass. Existing Fraction aggregation and normal-centering checks remain.
All RNA predictions, reference scales and decisions lock before DNA scoring.

The same48 tumour patients,20 normal rows and24/24 historical patient split
remain. Correlation is informative in23 comparison patients; loss recall in13.
All27 preceding method scores reproduce exactly, with30 methods total.
No patient disappears because a method fails or emits no classified calls.

| Method | Pearson(23) | Loss recall(13) | Classified precision | Loss FPR(23) |
|---|---:|---:|---:|---:|
| Gaussian regions,uniform median |.744107|.463553|.263445(17)|.044992|
| Gaussian regions,weighted median |.760619|.540476|.384804(17)|.042373|
| Gaussian regions,weighted mean |.755013|.540476|.345588(17)|.052557|
| Prior centered weighted median |.756351|.577106|.432456(16)|.042044|
| Installed native median |.738056|.495375|.348000(15)|.037392|
| Strong pseudobulk / own partition |.730261|.399679|.358571(14)|.041244|

Weighted median has20 wins/3 losses against strongest pseudobulk,18/5 against
native, and9 wins/12 losses/2 ties against the prior centered median. A few
larger gains yield the positive mean despite fewer wins against that preceding
candidate. Uniform/mean fail the declared positive-point prerequisite against
the preceding weighted median; retain them as negative/matched arms. Classified
precision has different call populations and is not unconditional accuracy.
[All primary results](../dev/reference_gaussian_regions_20261109/bulk_v1/report.json).

## Reference uncertainty rejects a promotion claim

The [uncertainty design](../dev/reference_gaussian_regions_20261109/UNCERTAINTY_DESIGN.md)
precedes resampling. Keep all three candidate arms, including those failing
point continuation. Family42 comprises36 candidate-versus-prior-baseline
contrasts and six weighted-versus-new-uniform contrasts. Endpoints are Pearson,
loss recall and loss FPR. Prior baselines are centered weighted median,native,
strongest pseudobulk and HD precision path mean. Matched uniform separates
weighting from changing boundaries; the preceding median isolates boundaries.

Original-reference000 reproduces EVERY primary candidate array and all30
patient score tables EXACTLY, gating eight resamples. All27 preceding methods
remain exact within each existing reference. Nested20,000 draws reuse the same
eight fixed reference draws and resample calibration/comparison patients with
the patient-balanced neutral-window threshold recalibrated each time. Seeds
2026110922/23.2583 exact threshold checks,6048 direct score checks, independent
fsum checks and Monte Carlo endpoint brackets pass. References do not multiply n.

| Weighted-median Pearson contrast | Paired gain |42-family interval|Holm p|Combined gate|
|---|---:|---|---:|---|
| Versus preceding centered median |+.004268|[-.013620,.046854]|1|Fails|
| Versus installed native |+.022563|[.007442,.075216]|1|Fails|
| Versus strongest pseudobulk |+.030358|[.007572,.078907]|1|Fails|
| Versus HD precision |+.027550|[.004163,.084291]|1|Fails|
| Versus matched new uniform |+.016512|[.007056,.062243]|1|Fails|

The rule still requires BOTH simultaneous lower bound>0 AND Holm p<.05.
Centred absolute-error bootstrap p-values and basic intervals are not dual
under these asymmetric reference-resampling distributions. Do not substitute
a positive interval for the failed test, shrink the family or claim added
biological replication. All42 results remain visible, including losses/FPR;
this development characterization does not control cumulative adaptive search
or validate repeated-workflow coverage. More draws cannot supply missing truth.
[Complete uncertainty](../dev/reference_gaussian_regions_20261109/uncertainty_v1/report.json).

## Compute and decision

| Jobs | Role | Account/QOS; CPUs / RAM / requested time | Outcome |
|---|---|---|---|
|23151063|independent numerics|nslab;1 /4GB /1h|complete0:0,24s|
|23151064|primary profile|innovation;1 /8GB /1h|complete0:0,52s|
|23151241|original-reference exact gate|innovation;1 /4GB /30min|complete0:0,52s|
|23151242–23151249|eight reference fits|alternating nslab/innovation;1 /4GB /30min|all complete0:0,53–58s|
|23151250|20,000-draw uncertainty|nslab;1 /8GB /1h|complete0:0,32s|

All execution uses CPU compute nodes006/009, a neutral directory and the
existing validated read-only noneditable wheel/native binary. Peak process
RSS251640KiB. Exact commands, versions, backend, dependencies, resources,
shared pre-edit baseline and immutable source_v1/v2 are retained in
[the task directory](../dev/reference_gaussian_regions_20261109/). Source_v1
runs primary/numerics; source_v2 adds the declared resampling mode. No package
test suite or GPU parity rerun is claimed for research-only helpers.

Reject a default replacement or strongest-baseline superiority claim. The point
profile improves but incremental boundary benefit is not established, and the
loss/FPR tradeoff does not uniformly improve. Former evaluations remain
regression-only. Suitable independent corresponding DNA, biological replication
and assay-specific truth resolution remain requirements. Clone recovery,
spatial validity and the5% whole-specimen false-split/power goal are unresolved.
