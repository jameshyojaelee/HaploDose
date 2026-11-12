# Native predictive RNA profile: performance candidate integration

The retained RNA-only smoother is now implemented as an opt-in native bulk
output. Default HMM/CNA/clone inference and uncertainty behavior are retained.
The promising prior181-line result remains historical regression evidence;
native integration alone does not establish a new accuracy gain.

## Change and scientific contract

`python/haplodose/smoothing.py::predictive_rna_profile` selects a running-median
width by conditional masked-gene absolute RNA prediction error. The estimand is
descriptive regional natural-log RNA abundance relative to the declared
reference. Feature selection, normalization and reference estimation are
conditioned on; they are not independent confirmation data. Reference mixture
and transcriptional-state effects can produce CNA-shaped RNA contrasts.
Independent corresponding DNA is needed for genomic interpretation.

`bulk.py::run_bulk(...,rna_profile_smoothing='predictive')` and CLI
`--rna-profile-smoothing predictive` add a float64 samples-by-bins array and
per-sample selection diagnostics after all existing inference has finished.
Units, axes, reference source and unvalidated status are explicit in run.json.
Defaults remain `none`. Arm/gap boundaries reset smoothing. Singleton genes
remain available for neighboring predictions; all-singleton inputs abstain
from selecting a width and retain their observed bin values.

The native gene-level input comes from declared raw counts. The exposed Kinker
prepared contrasts are transformed measurements and are never relabeled raw.
No calibrated uncertainty, absolute DNA CN or clone validation is inherited.

## Frozen checks and execution

[Design](../dev/native_predictive_profile_20261109/DESIGN.md) locks exact known
optima, independent scalar medians, singleton behavior, invalid-input checks,
native output parity and all181 retained-candidate identities. The new tests
are `tests/test_predictive_smoothing.py` and `tests/test_bulk_predictive.py`.
No existing test is weakened or changed by this implementation.

[Raw-bulk design](../dev/native_predictive_profile_20261109/RAW_BULK_DESIGN.md)
tests the original raw STAR panel, keeping48 tumour patients/20 reference
samples, the historical24/24 split, every12 comparator and common support.
RNA predictions freeze before scoring. Matched direct-gene widths1/5 isolate
smoothing from upstream changes. Loss recall, precision, achieved FPR, coverage
and correlation all remain visible; the same-cohort thresholds are descriptive.

The first pre-edit snapshot overlapped an edit and was rejected by its
consistency check (exit1, no SNAPSHOT_COMPLETE). It stays intact and unused.
The preceding verified source_v3 supplies the baseline; all75 existing native
and test files were checked, with only the intended bulk/CLI changes. New
snapshots source_v1/v2 completed before submission. The isolated environment
copies existing dependencies; shared installations and frozen outputs stay intact.

Environment job23136751 passed. Build/test23136752 built a wheel and passed4
Rust/411 Python tests but failed one new scalar-reference fixture: the private
rolling helper attempted NumPy indexing on a Python list. Source_v3 fixes this
by converting inputs to arrays; the fixture and every criterion are unchanged.
The failed build/test evidence and source remain intact. Four unstarted jobs
whose dependencies could never pass were canceled and replaced, not counted as
successful analyses. superseded_jobs.json records their pre-cancel states.

Corrected build23138003 passes412 Python tests and4 Rust tests;1 GPU-only test
is skipped on CPU. The installed selector reproduces every181 saved profile,
selected width and risk diagnostic exactly (23138004). This transports the
previous95 wins/33 losses/53 ties versus fixed-width-five pseudobulk into the
installed helper without treating the replay as new independent evidence.
Wheel SHA256 is `0e873bbaab5b98963255238f8b8d0914e443c4b0ce50fa0b109306dc6db8cf8e`;
native SHA256 is unchanged at
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.

## Raw-bulk performance: a useful correlation gain with a loss-recall cost

Job23138007 completes0:0. Every pre-existing native array/table matches the
frozen run exactly; all12 comparators and common support remain. The following
are historical comparison-split means, not prospective estimates:

| Method | Pearson (23 patients) | Loss recall (13) | Achieved loss FPR (23) |
|---|---:|---:|---:|
| Matched direct-gene median, width1 | .557408 | .304533 | .051795 |
| Matched direct-gene median, width5 | .630650 | .322711 | .045695 |
| Native predictive gene profile | .680245 | .311493 | .042193 |
| Naive smoothed | .598983 | .322482 | .040020 |
| Naive own partition | .730261 | .399679 | .041244 |
| Existing HD posterior | .720291 | .500183 | .057043 |
| HD path precision mean | .733069 | .577106 | .046896 |

Predictive versus matched width5 has15 correlation wins/4 losses/4 ties and
mean paired gain+.049595. Loss recall changes-.011218 (1 win/1 loss/11 ties).
Against naive smoothed, correlation improves+.081263 (19 wins/4 losses), but
recall changes-.010989. It trails naive own partition in correlation by.050016
and recall by.088187. The stronger existing HD path-precision comparator also
wins these endpoints. Therefore **retain the opt-in descriptive estimator;
reject a general bulk-accuracy or default-replacement claim**. No threshold is
retuned to repair this result.

All24 comparison patients remain in the output. Correlation/neutral truth is
informative in23 and loss truth in13. Precision remains conditional on making
classified calls: predictive has16 defined patient precisions, versus19 for
matched width5. Their separate means are.284772 and.200068; the actual paired
defined-case difference is+.047191 over16 patients. Do not interpret these
different call populations as an unconditional precision improvement. Every
method retains full classified callable coverage; undefined/no-call metrics,
all patient rows and counts are saved. [Full comparator results](../dev/native_predictive_profile_20261109/raw_bulk_v1/report.json).

Calibration thresholds use the same historical24 calibration patients. The
observed4.22% comparison FPR is a descriptive achieved rate, not a5% biological
guarantee or uncertainty coverage result. Shared references, uncertain DNA and
the exposed cohort remain limits. No new biological replication is supplied.

The native task has five executed CPU jobs: environment, failed first validation,
corrected validation, replay and raw bulk. Four complete0:0; the first validation
fails1:0 and stays recorded. The separate clone task has two completed jobs.
Four unstarted dependency jobs were superseded across both journals. Corrected checks run on AMD v6 nodes;
process/accounting resources and exact submission dependencies are preserved.
GPU parity was not rerun for this descriptive CPU export. The source/test
integrity record confirms existing scientific tests and clone code unchanged.
No frozen scientific result was overwritten. [Task records](../dev/native_predictive_profile_20261109).
