# Normal gene-bias correction improves point performance, not the pseudobulk verdict

A normal-centered reference-weighted regional median reaches mean development
Pearson.756351, versus.752184 preceding weighted median,.738056 installed native
median and.730261 strong pseudobulk. The36-contrast uncertainty family supports
an improvement over HD precision path means, but not over strongest pseudobulk
or installed native regions under the unchanged combined interval/Holm rule.
Keep the candidate experimental. No package default, original test, scientific
criterion or frozen output changes. This is an actual estimator experiment,
not another audit or an independent-validation claim.

## Scientific contract and smallest intervention

The estimand is a regional natural-log RNA effect relative to the typical normal
gene profile, conditional on the fitted programme expectation. The preceding
candidate uses the20 leave-one-out normal gene effects to estimate spread but
does not remove their gene-specific location. Pooled expectations and nuisance
fitting do not guarantee that each gene's median normal effect equals zero.

Compute each gene's median over those normal effects; subtract it from each
tumour's adjusted gene effect; recenter by the unchanged median of within-bin
gene medians. This is a stated reference-contrast change, not absolute DNA
dosage. Keep the same native RNA regions, genes, half-count transformation,
expected counts, genomic resets, normal-reference membership and heuristic
reference/count weights. Compare three fixed aggregators: uniform median,
reference-weighted median and reference-weighted mean. The uniform median is
the matched ablation separating location correction from reliability weighting.

Identification requires normal gene bias to transport to the tumour. Normal
tissue composition, transcriptional regulation, noisy offsets and reference
double-use can instead explain an apparent difference. Corresponding independent
DNA at common support tests profile prediction, not individual clone membership,
spatial assignment, absolute copy number, purity or ploidy. The same exposed
48 tumour patients,20 normal rows and historical24 calibration/24 comparison
split remain. Informative correlation n=23, loss-recall n=13. Reference resamples
reuse normal donors; they do not add biological replication. Former evaluations
remain regression-only. No verified unopened evaluation panel is introduced.

## Prespecified computation and preservation

[Primary design](../dev/normal_centered_regions_20261109/DESIGN.md) precedes
outcomes. Independent scalar normal-gene medians and two-stage centering,
reference-shift invariance and gene permutation pass, alongside the unchanged
Fraction weighted-median/mean and chromosome/arm/gap checks. Invalid/nonfinite
input contracts remain explicit. Gene effects and predictions are saved before
opening DNA scores. All24 preceding method scores reproduce exactly, including
the failed RNA-risk selector and harmful count-weighted mean. All27 methods,
patient rows, callable support and failure/undefined denominators stay visible.

Implementation extends the research
[reference helper](../benchmarks/hardening/reference_weighted_regions.py)
behind --normal-centering; default reference/count estimators retain their
previous calculations. No new production dependency, inference backend or
scientific threshold. Normal rows remain reference-only with missing candidate
estimates; they are not counted as tumour predictions.

| Method | Pearson(23) | Loss recall(13) | Classified precision | Loss FPR(23) |
|---|---:|---:|---:|---:|
| Normal-centered uniform median |.743489|.500183|.295139(16)|.046444|
| Normal-centered reference median |.756351|.577106|.432456(16)|.042044|
| Normal-centered reference mean |.751314|.577106|.438120(15)|.047033|
| Preceding reference median |.752184|.577106|.411623(16)|.042304|
| Installed native median |.738056|.495375|.348000(15)|.037392|
| Strong pseudobulk / own partition |.730261|.399679|.358571(14)|.041244|

Centered reference median has20 correlation wins/3 losses versus native,19/4
versus strong pseudobulk and16/7 versus preceding reference median. Its loss
recall matches the preceding weighted estimate, with somewhat lower FPR.
Against installed native, greater recall accompanies greater FPR. Classified
precision conditions on making calls; changing denominators cannot establish
an unconditional paired improvement. Full primary scores and all comparisons
remain in [bulk_v1](../dev/normal_centered_regions_20261109/bulk_v1/report.json).

## Conditional reference and patient uncertainty

The [uncertainty design](../dev/normal_centered_regions_20261109/UNCERTAINTY_DESIGN.md)
is frozen after point results and before resampling. All three candidates pass
the positive-point continuation rule against both native and strong pseudobulk;
that cheap gate is not promotion. The family36 contains three candidates x
four baselines x Pearson/loss recall/loss FPR. Baselines are installed native,
preceding reference median,strong pseudobulk and HD precision path mean.

Nine reconstruction jobs retain all24 preceding methods for each reference,
including rerunning the fixed failed RNA selector. Its new reference-refit
argument changes no default calculation. The
[refit driver](../benchmarks/hardening/normal_centered_refit.py) verifies that
reference000 reproduces EVERY primary candidate array and all27 patient score
tables exactly before dependent resamples. All eight existing normal-reference
resamples complete; they are sensitivity draws, not new patients or cohorts.

Reuse the exact patient-balanced neutral-window threshold and existing nested
bootstrap:20,000 draws, one of eight reference resamples per draw, independently
resampled calibration and comparison patients. Seeds2026110920/21. There are
2583 exact threshold checks and6048 direct score checks, plus independent fsum
effects and binomial Monte Carlo brackets on simultaneous endpoints. The
[uncertainty helper](../benchmarks/hardening/reference_weighted_uncertainty.py)
adds the declared36-family mode; preceding24-family modes keep their prior
seeds, methods and rules. This post-result family does not control the full
adaptive development history or establish repeated-workflow interval coverage.

| Centered reference median: Pearson contrast | Paired gain |36-family interval|Holm p|Combined criterion passes?|
|---|---:|---|---:|---|
| Versus installed native |+.018294|[.004673,.037792]|.059497|No|
| Versus preceding reference median |+.004166|[-.003783,.013980]|1|No|
| Versus strong pseudobulk |+.026089|[-.001757,.060047]|.475976|No|
| Versus HD precision path mean |+.023282|[.004229,.043865]|.010799|Yes, conditional development comparison|

The rule remains BOTH Holm p<.05 AND simultaneous lower bound>0. These summaries
are not dual under asymmetric bootstrap distributions. Do not use the favourable
native interval to override its failed Holm criterion, shrink the family, round
.059497 to significance or treat Monte Carlo refinement as new biological
evidence. The numerical brackets and finite-draw summaries stay visible.

Centered median loss-recall gain versus native is+.081731 with interval
[-.172887,.469256]; FPR improvement is-.004652 with[-.073313,.052406]. Versus
strong pseudobulk, recall gain+.177427 has[-.163225,.412196], and FPR improvement
-.000799 has[-.065952,.076631]. Neither loss claim passes. Centered uniform and
mean also fail superiority over strong pseudobulk; all36 contrasts are retained.
[Complete uncertainty](../dev/normal_centered_regions_20261109/uncertainty_v1/report.json).

## Compute and outcome

| Jobs | Role | Account/QOS; CPUs / RAM / time | Outcome |
|---|---|---|---|
|23150758|numerical fixtures|nslab;1 /4GB /1h|complete0:0,8s|
|23150759|primary estimator comparison|innovation;1 /8GB /1h|complete0:0,8s|
|23150854|original-reference exact gate|innovation;1 /4GB /30min|complete0:0,22s|
|23150855–23150862|eight reference reconstructions|alternating nslab/innovation;1 /4GB /30min|all complete0:0,12–15s|
|23150863|20,000-draw uncertainty|nslab;1 /8GB /1h|complete0:0,29s|

All12 jobs run on CPU compute nodes from a neutral directory using the same
read-only noneditable wheel/native binary and recorded imports/backend. Peak
process RSS249220KiB for the profile batch; uncertainty142320KiB. Scheduler
RSS0 is retained as reported, not substituted for process measurements. Exact
commands, dependencies, versions, resources and source hashes live in
[the task directory](../dev/normal_centered_regions_20261109/). Primary uses
the shared binomial task source_v1; resampling uses its own source_v2. Baseline
is the explicitly shared pre-edit binomial task snapshot. All snapshots and
historical reports remain intact. Full421 Python/4 Rust tests are not rerun for
research-only changes; no new full-suite/GPU claim is made.

The candidate provides a better point profile and a conditional advantage over
one strong HaploDose comparator. It does not yet establish pseudobulk superiority,
transportability or biological calibration. No promotion follows. Independent
corresponding specimens and suitable truth remain an evidence requirement that
additional computation cannot replace. Original clone recovery and the5%
whole-specimen false-split/power objective remain unresolved.
