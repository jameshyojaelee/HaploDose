# Bulk amplitude restriction on fixed RNA-selected boundaries

This experiment separates continuous regional RNA-effect estimation from the
HMM's discrete amplitude restriction. It does not change boundaries, nuisance
fits, references, callable support, calibration rules or old comparators. The
[frozen design](../dev/bulk_continuous_amplitude_20261109/DESIGN.md) was written
before execution. The completed independent audit supports a promising conditional
weighted-amplitude result, with no default promotion or superiority claim.

## Contract, methods and falsifiers

The estimand is a constant natural-log regional RNA effect relative to the
executed reference/programme expectation and fitted centre. On the existing
Viterbi runs, two Gaussian working estimators are ordinary and inverse-variance
means of `log((Y+.5)/(E+.5))-centre`; the latter uses fitted
`sd_j^2+1/(Y_j+.5)`. This is a conditional log-rate approximation, not exact PLN
maximum likelihood, calibrated confidence or a DNA copy-number estimate.

The assumption for genomic prediction is that regional RNA effect tracks DNA
after nuisance adjustment. Regional transcription/reference mismatch and removal
of useful shrinkage are competing explanations. Corresponding independent DNA
at matched resolution would distinguish predictive benefit; fixed-model fit
alone cannot establish it. The panel reuses48 raw-STAR COAD/READ patients and20
reference donors with24/24 exposed roles. ABSOLUTE GRCh37 truth remains uncertain
and gene anchored to hg38 RNA, not lifted over. There is no new biological
replication or independent evaluation.

All RNA arrays are completed before scoring DNA is loaded. All48 original
posterior/effect/path fits must reproduce exactly. All10 historical methods are
retained, including strong naive own-partition, and the two new estimators are
scored on identical windows. Patient-equal neutral5% calibration is retained,
with achieved comparison FPR, recall, precision and correlation reported together.
This is different from the whole-specimen clone false-split target.

Independent checks include scalar integer count aggregation, partition coverage
and resets, analytic/root/Decimal80 amplitudes, translation/precision-rescaling
invariants, gene-weighted window medians, exact Fraction calibration ranks and
patient-level endpoint recounts. Reference recalibration checks numerical call
stability; unresolved boundaries stay unresolved. Paired effects retain matched
patient denominators. No new hypothesis tests or superiority claim are planned.

No improvement on fixed boundaries rejects amplitude restriction alone as a
sufficient repair. Higher recall with worse error or precision is not success.
No default is promoted from exposed point estimates. A promising result would
need paired uncertainty, reference robustness and independent frozen evaluation.

## Failure and execution journal

| Job | Role | Recorded outcome |
|---|---|---|
|23124626|First48-patient panel|Failed42s,exit1:0 after fitting/scoring,at final report:missing `write_json` import. Outputs retained,not accepted as successful run.|
|23124627|First audit|Cancelled while pending,before execution,after audit syntax error was caught.|
|23124641|Syntax-corrected audit|Cancelled while pending after panel dependency failed; no analysis executed.|
|23124667|Corrected panel into `panel_v2`|Completed16s,exit0:0; all48 original fits match exactly.|
|23124668|Independent audit into `audit_v3`|Failed15s at window precision comparison; scalar bin/segment checks had passed but the reference omitted the scorer's float32 output cast.|
|23124691|Storage-aware independent audit into `audit_v4`|Completed16s,exit0:0; exact cast reproduction plus unchanged amplitude checks and refined decisions.|

The initial static audit parse caught an unmatched parenthesis in a Decimal list;
the submission command incorrectly continued after that failure. The panel
script itself had parsed. The pending audit was cancelled and source v2 corrected
the syntax. Source v3 also adds the missing panel import and an explicit audit
panel-directory selector. Subsequent submission requires successful AST parsing
before proceeding. These are execution repairs, not scientific criterion changes.
All sources, failed outputs and exact commands remain under
[bulk_continuous_amplitude_20261109](../dev/bulk_continuous_amplitude_20261109).

The corrected panel requests innovation CPU4/8GB/2h; audit requests nslab
CPU1/8GB/1h. Both use the unchanged403-Python/4-Rust validated noneditable wheel.
No production file, original test, dependency, scientific threshold, shared
installation or frozen scientific result is edited by this experiment.

The storage audit failure is retained. The existing `subarm.window_median` writes
float32 medians; the first reference compared its float64 values against that
storage using only the tighter amplitude arithmetic guard. The final audit
retains that amplitude guard, independently reproduces the float32 arrays exactly,
and separately accounts for the declared half-ULP storage rounding. It also
recalibrates the full-precision predictions and checks every callable decision.
No scorer cast, threshold, scientific criterion or prediction was changed to pass.

## Independent checks and descriptive outcomes

Audit v4 verifies54288 scalar integer bin sums,2056 complete run partitions,
96 Decimal80 first/last-region checks, scalar score-root amplitude estimates,
translation/weight-rescaling invariants and2304 patient endpoints across12 methods.
All10 old method arrays/scores, truth, support and sample order remain identical.
The maximum recorded window discrepancies are5.69e-8/5.65e-8 for ordinary/weighted
amplitudes, entirely consistent with the explicit float32 storage step. The
maximum segment arithmetic guard is3.94e-12, unchanged from the failed audit.
Full-precision reference recalibration changes0 callable decisions for both arms.

There are24 comparison patients throughout. Correlation is defined for23, loss
recall for13 and neutral FPR for23. Precision is undefined when a method makes no
classified loss calls: it is defined for14 HMM patients,16 ordinary-amplitude
patients and15 weighted-amplitude patients. These different subsets must not be
compared as though their marginal precision means were a paired effect.

| Exposed comparison method | Mean Pearson,23 patients | Mean loss recall,13 patients | Mean loss FPR,23 patients | Mean precision,method-specific callers |
|---|---:|---:|---:|---:|
| HMM per-bin posterior |0.720291|0.500183|0.057043|0.411818 (14)|
| Same HMM partition,ordinary mean |0.722760|0.572299|0.049631|0.416842 (16)|
| Same HMM partition,precision-weighted mean |0.733069|0.577106|0.046896|0.458667 (15)|
| Strong naive own-partition |0.730261|0.399679|0.041244|0.358571 (14)|

The paired effects below use exactly the patients evaluable for both arms in
each endpoint. Positive correlation/recall/precision and negative FPR favor the
new method. Precision remains conditional on both methods making classified calls;
this does not establish full-cohort precision benefit or erase abstention.

| New estimator | Comparator | Pearson,n23 | Recall,n13 | FPR,n23 | Paired precision,n |
|---|---|---:|---:|---:|---:|
| Ordinary mean | HMM posterior |+0.002469|+0.072115|−0.007411|−0.006856,n14|
| Weighted mean | HMM posterior |+0.012778|+0.076923|−0.010147|+0.008182,n14|
| Ordinary mean | Naive own-partition |−0.007502|+0.172619|+0.008387|−0.001323,n13|
| Weighted mean | Naive own-partition |+0.002807|+0.177427|+0.005652|+0.014872,n13|

All paired IDs, missing denominators and conditional interpretations are retained
in `audit_v4/paired_effects.json`; all patient call/coverage counts remain in
`panel_v2/patient_scores.tsv`. The apparent ordinary-mean precision gain in the
marginal table reverses on paired patients. Do not claim a precision improvement
for it. Earlier descriptions of24 comparison patients refer to total population,
not a guarantee that24 contribute to every endpoint.

[Observed pattern, not a guarantee] Precision weighting improves the old HMM's
four paired point endpoints on this fixed-reference exposed comparison. It has
slightly higher correlation and higher recall than the strong naive baseline,
but also higher achieved FPR. A favorable mean FPR below.05 is not a confidence
bound, prospective calibration or the5% specimen clone guarantee. No uncertainty
interval or significance test is estimated in this bounded mechanism study.

## Decision and next evidence gate

The fixed-boundary result supports investigating continuous RNA-amplitude
estimation further. It does not establish a globally better segmentation method,
DNA-state inference, purity or general package superiority. The improvement is
conditional on existing HMM-derived boundaries and nuisance/reference estimates.
The simpler ordinary mean does not show the same paired benefit.

Retain both arms as research-only evidence. Before any promotion, freeze a
reference-robustness and biological-unit uncertainty study using the existing
reference resamples and complete fixed comparator family. Do not tune centering,
weights or boundaries on the exposed outcomes. Only an independently frozen
suitable evaluation can establish a broader performance claim. Clone recovery,
false splitting, membership calibration and spatial validity remain separate,
unresolved requirements; this bulk experiment does not count as their repair.

All six job records remain:2 completed,2 failed,2 cancelled before execution.
Source snapshots v1–v4 and their unchanged inventories remain alongside the failed
outputs. `accounting_final.tsv`, `source_integrity_checkpoint.json`, `jobs.json`
and `analysis_JOB` preserve status, hashes, commands, backend/import identity and
process resources. No full production suite was rerun for research-only code;
the validated production artifact remains403 Python/4 Rust passes.
Successful panel/audit process peaks were946380/188928KiB, with elapsed15.17/14.55s;
the scheduler records16s each. No allocation exceeded8GB. Source inventories
retain744,744,744 and745 unchanged file hashes for versions1–4, respectively.
