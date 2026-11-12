# Performance candidates: end-to-end recovery and continuous profiles

Two scientific candidates were implemented and run on compute nodes in response
to the request to prioritize performance. Neither establishes the requested win
over the strongest baseline or fixes the three recovery regressions. Production,
frozen results, original recovery assertions and false-split criteria are unchanged.
These are performance experiments, not new biological validation.

## End-to-end clone candidate

[Implementation](../benchmarks/hardening/performance_candidates.py) replaces
the working difference model with the previously fixed continuous regional-effect
model throughout actual proposal confirmation, merging and export calls. It
integrates a signed normal RNA effect shared within each run, centered on the
existing step with SD equal to that step, and retains the same state-path prior.
Proposal-half assignments remain; no full-count membership substitution occurs.
Every-bin .95 support, core lengths, group support and all six recovery criteria
remain unchanged. Estimand is conditional assigned-pool RNA contrast, not a DNA
copy-number difference. Regional transcriptional state remains an alternative;
same-cell DNA would distinguish genomic identity.

The [frozen design](../dev/performance_candidates_20261109/DESIGN.md) precedes
outcome inspection. Job23130970 executes all five original full workflows and
1435 testable difference comparisons in6m52s. Recovery remains **2/5**:

- 20261031 and20261107 pass all original criteria.
- 20261108 still returns one S1 tumour group, assignment accuracy.5960 and normal
  FPR.1111, failing the unchanged thresholds.
- 20261109 still returns one S1 group, accuracy.6020, no chr12 distinction.
- 20261110 still returns one S1 group, accuracy.5996, no chr12 distinction.

[Full report and all criteria](../dev/performance_candidates_20261109/recovery_v1/report.json),
[candidate calls](../dev/performance_candidates_20261109/recovery_v1/model_calls.json)
and original public output tables remain available. A slight accuracy change in
an already passing case is not a recovery milestone. No large null/power batch
was launched for this unsuccessful recovery candidate; no promotion follows.

### Preceding actual-merge performance probe

Before the pause, jobs23130857/23130858 captured all4 actual merge contrasts and
tested the same prior. Exact old pool/public-output/30-criterion parity passes.
Independent point-limit HMM and complete-chr12 quadrature errors are at most
3.31e-14 for the quadrature check. In20261108, one failed merge reaches only
.942313 support and no core; the other gains a qualifying209–216 core with
minimum support.993446. This partial conditional result does not predict full
workflow recovery, as the subsequent end-to-end experiment demonstrates.

[Merge-stage design](../dev/merge_amplitude_20261109/DESIGN.md),
[capture report](../dev/merge_amplitude_20261109/capture_v1/report.json),
[all four conditional results](../dev/merge_amplitude_20261109/panel_v1/comparisons.json)
and [source](../benchmarks/hardening/merge_amplitude.py) are retained.
The .9423 result was not accepted by lowering .95 or selecting the favorable merge.

## Continuous regional profile candidate

The new estimator removes fixed copy-number state amplitudes: exact weighted
optimal partition fits piecewise constant RNA effects to existing RNA tracks,
respecting genomic resets. It minimizes weighted residual SSE plus
3*log(n)*(number_of_segments-1) per block. This is an MDL/BIC regularization
approximation: log(n) for a new continuous mean and2*log(n) for an unspecified
breakpoint. It is not a calibrated CNA discovery threshold. The Gaussian variance
and reference are working assumptions; correlated regulation remains a biological
alternative. Corresponding DNA at matched resolution tests prediction accuracy.

No DNA-selected penalty grid was used. Predictions precede scoring; genes,
callable windows, reference construction,24 calibration and24 comparison patients
and the full12-comparator family are fixed. The panel is exposed development
evidence, not prospective evaluation. One comparison patient lacks usable
correlation/neutral truth, and only13 have loss truth. Neither patients, windows
nor reference resamples are interchangeable biological replication units.

Job23130971 fits the first reference;23130997 applies the identical candidate to
all9 existing reference constructions. Job23131131 tests coordinatewise median
reference prediction for **every** method, including the strong naive baseline.
[Reference follow-up](../dev/performance_candidates_20261109/REFERENCE_FOLLOWUP.md)
and [ensemble design](../dev/performance_candidates_20261109/ENSEMBLE_DESIGN.md)
record when these follow-ups were defined. No reference subset/weight tuning.

| Comparison-panel marginal mean | Candidate | Strong naive own-partition |
|---|---:|---:|
| Single-reference Pearson |.731090|.730261|
| Single-reference loss recall |.613507|.399679|
| Single-reference achieved loss FPR |.044244|.041244|
| Median-ensemble Pearson |.732073|.731508|
| Median-ensemble loss recall |.618315|.399679|
| Median-ensemble achieved loss FPR |.047800|.041244|

These are marginal descriptive means, not confidence bounds or independent wins.
Across references, recall improves9/9, correlation improves2/9 and FPR is
higher8/9. The candidate does not robustly beat the strong naive correlation.
The precision-weighted HMM-path estimator also has higher correlation than this
candidate in the original and ensemble comparisons. The method is not promoted.

## Paired performance assessment

Job23131272 characterizes two candidate arms ×12 comparators ×4 endpoints=96
contrasts. Unit is patient. Exact one-sided sign tests target the median paired
improvement, with Holm over96; simultaneous one-sided95% lower median bounds use
exact binomial order statistics and Bonferroni96. Zeros count as nonpositive.
Mean paired effects are reported separately; the sign test does not test means.
Missing paired endpoints and excluded patients remain explicit. Undefined
precision is not silently treated as zero. The
[post-hoc assessment design](../dev/performance_candidates_20261109/EFFECTS_DESIGN.md)
does not relabel this development-selected analysis as independent confirmation.

Against strong naive, paired mean recall improves.213828(single) and.218636
(ensemble), but median paired recall change is0 on13 patients. Neither candidate
has a median superiority decision against strong naive on any endpoint. The
paired correlation mean differences are only.000829/.000565, with simultaneous
lower median bounds-.062891/-.064863. Seven descriptive median superiority
decisions concern correlations versus weaker hd_lr, naive or smoothed exports.
They do not satisfy the strongest-comparator objective.

Marginal precision means initially look higher, but defined-pair mean precision
differences against strong naive are-.000277(single) and-.000114(ensemble), on13
paired patients; median changes are0. Different undefined-precision populations
explain why marginal means cannot be used as paired improvements.

[All96 paired effects](../dev/performance_candidates_20261109/effects_v1/effects.tsv),
[all patient scores](../dev/performance_candidates_20261109/profile_v1/patient_scores.tsv),
[reference results](../dev/performance_candidates_20261109/profile_sensitivity_v1),
[ensemble results](../dev/performance_candidates_20261109/ensemble_v1/summary.tsv)
and [effect implementation](../benchmarks/hardening/performance_effects.py)
retain all outcomes. A population mean-effect uncertainty claim was not established.
No Kinker pseudobulk win, spatial win or biological calibration follows.

## Execution and decision

Five performance jobs complete0:0: recovery6m52s, initial profile16s, nine-reference
screen27s, fair ensemble10s and paired effects3s. Earlier conditional merge jobs
complete0:0 in1m14s/14s. Actual process memory and scheduler values are retained
separately under the two task directories. Source_v1 executes recovery/initial
profiles; source_v3 executes reference sensitivity; source_v4 the ensemble;
source_v5 paired effects. Source_v2 is an unused intermediate. All use the
validated noneditable wheel, isolated environment and immutable source snapshots.

[Performance commands](../dev/performance_candidates_20261109/jobs.json),
[accounting](../dev/performance_candidates_20261109/accounting.tsv),
[resources](../dev/performance_candidates_20261109/resource_summary.json),
[conditional merge commands](../dev/merge_amplitude_20261109/jobs.json) and
[conditional merge resources](../dev/merge_amplitude_20261109/resource_summary.json)
are preserved. CPU resources only; no arrays, interactive QOS or shared installs.

Keep these implementations research-only. Neither is a justified default or an
accomplished main performance milestone. The next recovery candidate must address
the observation/membership information rather than reuse full-count selected
contrasts. Profile work must improve paired outcomes versus the strong baseline
under reference variation, with uncertainty and achieved error assessed together.
Frozen historical results and all scientific criteria remain intact.
