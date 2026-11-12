# Continuous RNA amplitudes under donor-reference resampling

This follows the [fixed-reference amplitude study](BULK_CONTINUOUS_AMPLITUDE_2026-11-09.md)
without changing its two estimators, ten comparators, selected boundaries or
scientific criteria. The [prespecified design](../dev/bulk_amplitude_reference_20261109/DESIGN.md)
requires exact original-reference reproduction before running the eight existing
normal-donor resamples. These are sensitivity analyses on exposed patients,
not new biological cohorts or independent evaluation.

## Scientific contract and decisive comparison

The estimand is the regional natural-log RNA effect conditional on the fitted
reference/programme expectation and centre, estimated by an ordinary or
inverse-variance mean inside existing HMM Viterbi regions. Its observation model
is the same Gaussian log-rate working approximation used in the preceding
experiment. The variance is fitted bin variance plus count-sampling variance.
Neither amplitude is an absolute DNA dosage or a calibrated interval.

The key assumption is that the reference/nuisance fit leaves a genomic dosage
contrast distinguishable from tissue composition and regional transcription.
The strongest alternative to an apparent predictive gain is reference-dependent
amplitude/shrinkage behaviour on these already exposed patients. The smallest
available check holds all other choices fixed while replacing only the frozen
donor-bootstrap reference. Selecting the reference with the best DNA score is
prohibited. Independent corresponding DNA and suitable biological replication
remain necessary to validate the genomic interpretation at the declared
approximately 10Mb gene-anchored resolution; assembly correspondence remains a
truth limitation, not something reference resampling resolves.

All 48 tumour patients, 20 reference rows, 24 calibration patients and 24 exposed
comparison patients remain fixed. The eight bootstrap references reuse those
normal donors. They do not multiply biological n. All 12 methods are retained:
the original five bulk exports/baselines, posterior effect, ideal constant,
physical-distance and shuffled-distance HMM controls, naive own-partition,
and the two continuous amplitudes. Historical naive/naive-smooth log2 units
remain historical; thresholds are calibrated separately for each score.

For each reference, independently recount all patient endpoints, then compare
each candidate with every historical comparator for Pearson correlation, loss
recall, neutral-bin FPR and classified-call precision. This is 80 descriptive
paired contrasts per reference. Lower FPR is favourable; the other three
endpoints favour larger values. Precision conditions on both methods making
classified calls, so its paired population may change with the reference.
Undefined endpoints and changing denominators remain explicit.

Reversal of a favourable point effect under reference resampling falsifies
reference-invariance of that effect. Uniformly favourable point effects alone
do not establish statistical superiority, confidence coverage or biological
validity. Neither estimator is promoted from this experiment. Paired patient
uncertainty with recalibrated thresholds remains a separate required gate.

## Computational gates and implementation

[Primary runner](../benchmarks/hardening/bulk_amplitude_reference.py) reconstructs
each frozen fit with the validated wheel. Every tumour must exactly reproduce
the old float32 log-ratio export, fitted centre/beta/tau, segment boundaries,
states and segment posterior mean effects. Integer bin totals are conserved.
The original-reference run must additionally reproduce all 12 preceding window
arrays and patient scores exactly. All RNA predictions finish before DNA
scoring values are loaded. Truth, callable support and sample identities must
remain identical across references.

[Independent audit](../benchmarks/hardening/bulk_continuous_amplitude_audit.py)
retains the existing scalar integer aggregation, long-double log-rate/variance,
SciPy scalar score-root, Decimal80, reset/partition, sorted gene-window median,
exact Fraction threshold and patient endpoint checks. The reference-mode option
also verifies every tumour column and every resampled normal column against its
recorded original specimen, including donor/role correspondence. It retains the
binary64 arithmetic guards and explicitly checks unchanged float32 window
storage. Full-precision recalibration records any callable decision differences.

The geometry controls retain the independently studied maximum posterior
discrepancy criterion of 1e-4 against float64 recursion. There is no new
universal epsilon or relaxed test. The
[summary](../benchmarks/hardening/bulk_amplitude_reference_summary.py) requires
all nine successful independent audits and verifies their paired effects from
the patient tables before reporting every reference, comparator and denominator.

## Execution and provenance

The [job journal](../dev/bulk_amplitude_reference_20261109/jobs.json) records all
19 submitted jobs: original-reference gate 23124870; eight resample fits
23124911–23124918; nine dependent audits 23124919–23124927; final dependent
summary 23124928. No arrays, interactive QOS, installations or production edits.
Fits use innovation/innovation, cpu, v6&amd, 4 CPUs, 8GB, 2 hours; audits and
summary use nslab/nslab, cpu, v6&amd, 1 CPU, 8GB, 1 hour. Each job records its
actual backend, isolated imports, wheel/native/source identity, environment,
runtime, resource usage and exit status.

The original-reference gate uses source_v1; subsequent jobs use source_v3.
source_v2 was retained but never executed: static inspection found that the
original-reference sample table lacks the bootstrap-only source_sample field;
the audit now explicitly falls back to the original sample identifier. No
failed scientific check was bypassed. The primary runner is unchanged across
these snapshots. Production wheel/native identity remains the
[exact-profile artifact](../dev/exact_profile_merge_20261109/validation_23121461/installed_artifact.json).
Frozen results are read only; new artifacts live under the new task directory.

## Results and decision

All 19 jobs completed with exit 0:0. All nine references reproduce their frozen
fits; the original reference additionally reproduces all preceding 12-method
outputs and scores exactly. The independent audits total 488,592 scalar bin
sums, 18,649 regions, 864 Decimal80 checks, 20,736 patient endpoints and 180
normal-column source correspondences. Full-precision recalibration changes zero
callable decisions across both estimators and all nine references. These counts
are computational checks, not independent biological sample sizes.

The [full 80-contrast reference table](../dev/bulk_amplitude_reference_20261109/summary_v1/reference_stability.tsv)
retains every result. The table below shows paired weighted-amplitude effects;
ranges are the minimum and maximum across the eight reference draws, not
confidence intervals. Favourable counts use lower FPR and higher other endpoints.

| Comparator | Endpoint | Original effect | Reference-draw range | Favourable draws | Paired patients per draw |
|---|---|---:|---:|---:|---:|
| HMM posterior | Pearson | +0.012778 | +0.009168 to +0.020879 | 8/8 | 23 |
| HMM posterior | Loss recall | +0.076923 | -0.004808 to +0.076923 | 5/8; one tie | 13 |
| HMM posterior | Neutral FPR | -0.010147 | -0.013846 to +0.011158 | 5/8 | 23 |
| HMM posterior | Precision | +0.008182 | -0.011921 to +0.060637 | 6/8 | 13–15 |
| Naive own-partition | Pearson | +0.002807 | -0.006675 to +0.003772 | 5/8 | 23 |
| Naive own-partition | Loss recall | +0.177427 | +0.084707 to +0.177427 | 8/8 | 13 |
| Naive own-partition | Neutral FPR | +0.005652 | -0.008289 to +0.021864 | 1/8 | 23 |
| Naive own-partition | Precision | +0.014872 | -0.041648 to +0.068469 | 5/8 | 13–14 |

[Observed pattern, not a guarantee] Weighted amplitudes retain a correlation
advantage over posterior amplitudes in this finite reference bank, but the
other apparent HMM improvements reverse in some draws. Against strong naive
segmentation, higher recall persists while FPR is worse in seven of eight
draws; the small correlation advantage reverses in three. Ordinary amplitudes
have lower correlation than naive own-partition in all eight draws, despite
higher recall in all eight. Changing precision populations cannot support an
unconditional full-cohort precision claim.

This rejects the interpretation that the earlier fixed-reference improvement
was a stable joint improvement in correlation, recall, FPR and precision. It
does not show that continuous amplitudes have no useful role. The next
uncertainty analysis must resample calibration patients, comparison patients and
the reference bank jointly, retain all 80 contrasts, and keep undefined/caller
populations visible. It must not treat these eight draws as independent cohorts
or choose favourable references. No confidence interval, new hypothesis test,
independent-validation claim or default promotion has been made.

Resource records distinguish scheduler sampling from process high-water marks.
Fits took 33.81–53.82 seconds of process wall time and 960,416–980,264 KiB peak
process RSS; audits took 12.94–15.04 seconds and 195,632–203,684 KiB. The summary
took 2.45 seconds and 82,208 KiB. Some short jobs have zero scheduler MaxRSS;
their `/usr/bin/time` records supply the process measurement. See
[accounting](../dev/bulk_amplitude_reference_20261109/accounting_final.tsv),
[resources](../dev/bulk_amplitude_reference_20261109/resource_summary_final.json),
[audit inventory](../dev/bulk_amplitude_reference_20261109/audit_inventory.json)
and [source integrity](../dev/bulk_amplitude_reference_20261109/source_integrity_checkpoint.json).
All 75 production/native/test source files match the immutable snapshot. The
full production test suite was not rerun for this research-script-only change.
