# Joint reference, calibration and patient uncertainty for bulk amplitudes

The [reference-sensitivity study](BULK_AMPLITUDE_REFERENCE_2026-11-09.md) found
that several fixed-reference gains reversed when the normal reference changed.
This stage estimates uncertainty without choosing a favourable reference,
candidate, comparator or endpoint. It retains both continuous-amplitude
estimators, all ten historical controls and all four endpoints: 80 contrasts.
The [design](../dev/bulk_amplitude_uncertainty_20261109/DESIGN.md) was frozen
before the pilot and before any uncertainty outcome was viewed.

## Estimand, assumptions and limits

Each statistic is a paired comparison-patient mean endpoint difference.
Pearson, loss recall, neutral-bin FPR and classified-call precision have
different eligible populations. Precision conditions on both methods making
classified calls; it is not unconditional cohort precision. Every draw retains
the paired population size. Missing endpoints remain undefined, not zero.
Effect orientation is positive for improvement, including a reversed sign for
FPR; the final report also retains the original signed point effect.

The working resampling model draws uniformly from the eight existing
normal-donor bootstrap fits, resamples 24 calibration patients, recalculates
the exact equal-patient 5% neutral threshold for every method, then resamples
24 comparison patients. Reference 000 supplies the point estimates. The
20 normal-reference donors and 48 tumour patients have no patient-ID overlap;
the 24/24 tumour split contains no repeated patient. Independence of the two
donor populations is an assumption supported by these identifiers, not proof
that unrecorded laboratory or sampling dependence is absent.

The eight references are a finite bank of resamples of one donor pool, not
independent cohorts. Intervals condition on that bank and the exposed RNA/DNA
measurement chain. Uncertain assembly correspondence and old gene-anchored DNA
truth remain limitations. Regional transcription and reference composition
remain alternatives to genomic dosage. A precise bootstrap cannot distinguish
observationally equivalent biological explanations or establish empirical
coverage. Independent corresponding DNA and suitable biological replication
are still needed at the claimed resolution.

## Exact computation before scaling

The [runner](../benchmarks/hardening/bulk_amplitude_uncertainty.py) introduces
research-only caches for the existing threshold and endpoint definitions.
Neutral calibration uses arbitrary-precision integer cumulative masses with a
common denominator, so the equal-patient quantile remains exactly 1/20.
Multiplicity of resampled patients changes their mass explicitly. Empty
neutral contributors are excluded in the same way as the existing exact
scorer. Strict score-threshold comparisons and finite classified-call support
remain unchanged. No approximation, threshold jitter or convenient epsilon is
introduced for speed.

The [independent audit](../benchmarks/hardening/bulk_amplitude_uncertainty_audit.py)
compares calibration with Fraction arithmetic and endpoints with direct scalar
truth counting. Fixtures include unequal support, ties, duplicated patients,
empty support, nonfinite scores and thresholds. All indexed resampling
identities are replayed. Selected actual draws are independently reconstructed
with compensated sums and an operation-count rounding bound. Original point
effects are anchored to the preceding independent 80-contrast audit.

Pilot job 23125800 and audit job 23125801 both completed with exit 0:0.
The audit passed 274 fixture checks, all 1,000 indexed draw identities, 780
Fraction thresholds, 112,320 scalar patient endpoints, 5,200 scalar paired
effects and all 80 original point effects. The pilot is retained separately
and is not counted twice in the million-draw result. Pilot process wall time
was 11.07 seconds with 1,834,704 KiB peak RSS; audit time was 17.53 seconds with
1,833,008 KiB. This is computational validation, not confidence coverage.

## Fixed numerical and scientific gates

The fixed budget is one million draws, RNG default_rng(20261109), split across
20 independent job submissions with disjoint contiguous 50,000-draw ranges.
Every chunk follows the same indexed stream and saves all 80 effects, all paired
population sizes, 12 thresholds and the reference/calibration/comparison IDs.
The independent audit checks every chunk before the
[summary](../benchmarks/hardening/bulk_amplitude_uncertainty_summary.py) can run.
Every endpoint and failed/untestable contrast remains in the declared family.

Centered approximate bootstrap p-values and basic intervals refer to the same
paired mean effect. Holm uses the full 80-member family; simultaneous basic
intervals use Bonferroni with alpha=.05. Constant point-pair or bootstrap
distributions are untestable. Any undefined bootstrap statistic prevents an
unconditional gate for that endpoint; conditional defined-draw summaries and
the undefined rate remain visible. No arbitrary valid-draw fraction is used
to hide failures.

The family-interval tail is .05/(2*80)=.0003125, giving approximately 312.5
draws per extreme tail at the fixed budget. That motivates the numerical
budget, not its scientific adequacy. Simultaneous numerical bounds cover
80 p-values and 160 family-interval endpoints, a family of 240 quantities.
They use exact binomial/Clopper-Pearson bounds and binomial order-statistic
quantile bounds with .05/(2*240) per numerical tail. Independent rational
quantile interpolation and a separate Holm implementation check the report.
If a numerical bound crosses its decision boundary, the decision remains
unresolved. Passing a conditional bootstrap gate does not validate biological
coverage or authorize a production default.

## Execution record and current status

The [journal](../dev/bulk_amplitude_uncertainty_20261109/jobs.json) records every
submission. Source_v1 contains the pilot and audit. Source_v2 adds the final
summary; pilot and full-run calculation/audit source bytes are identical.
All jobs use the noneditable validated wheel in its isolated environment and
record source/import/native identity, backend, runtime, RSS and failures.
No production dependency, inference code, frozen result or scientific
criterion is changed.

The 20 draw jobs are 23126257,23126259 through 23126295 in steps of two;
the corresponding dependent audits are 23126258 through 23126296 in steps
of two. Summary 23126297 depends on all 20 audits. Each requests cpu,
v6&amd, 1 CPU, 8GB, 2 hours, alternating innovation/nslab for draws and the
opposite account for their audits. No arrays or interactive QOS.

All 44 submitted jobs completed with exit 0:0: pilot/audit, 20 chunks, 20
chunk audits and two summary versions. The accepted report is
[summary_v2](../dev/bulk_amplitude_uncertainty_20261109/summary_v2/summary.json),
job 23126610. The first summary remains preserved. Static review found that its
input-record names for original-reference audits collided with names for the
new chunk audits. The corrected collector prefixes reference-bank records
separately. Both scientific JSON reports are byte-identical; all 20 accepted
chunk-audit paths are now independently checked. This was a provenance repair,
not a change to data, effects, tests or numerical criteria. See
[provenance check](../dev/bulk_amplitude_uncertainty_20261109/summary_provenance_check.json).

The 20 full audits verify every one of the million resampling identities,
5,480 fixture checks, 15,600 exact Fraction thresholds, 2,246,400 scalar patient
endpoints, 103,980 defined scalar paired effects and 1,600 original point
effect anchors. Undefined selected-draw effects are explicitly checked as
undefined, not silently counted as passing numerical comparisons. These are
computational checks, not independent biological replicates.

## Scientific result and decision

Six of the 80 conditional contrasts clear the declared Holm/family-interval
gate even under simultaneous numerical bounds. They are Pearson correlation
for each candidate versus hd_lr, naive and naive_smooth. Neither candidate
establishes superiority over posterior HMM profiles or strong naive
own-partition segmentation. The negative decision is lack of established
improvement, not proof of equivalence or proof that amplitude estimation has
no useful role.

For the weighted candidate, the table reports effects oriented so positive
favours the candidate, including the reversed FPR sign. Intervals are the
approximate simultaneous 80-family basic-bootstrap intervals, not validated
biological confidence bounds. They condition on the finite reference bank.

| Comparator | Endpoint | Oriented point effect | 80-family interval | Gate |
|---|---|---:|---:|---|
| HMM posterior | Pearson | +0.012778 | [-0.020067, +0.043301] | Not established |
| HMM posterior | Loss recall | +0.076923 | [-0.431161, +0.382800] | Not established |
| HMM posterior | Neutral FPR | +0.010147 | [-0.055722, +0.076713] | Not established |
| HMM posterior | Precision | +0.008182 | Conditional defined-draw summary only | Untestable: 708 undefined draws |
| Naive own-partition | Pearson | +0.002807 | [-0.032275, +0.037715] | Not established |
| Naive own-partition | Loss recall | +0.177427 | [-0.179879, +0.588187] | Not established |
| Naive own-partition | Neutral FPR | -0.005652 | [-0.076845, +0.097410] | Not established |
| Naive own-partition | Precision | +0.014872 | Conditional defined-draw summary only | Untestable: 121 undefined draws |

Sixteen contrasts have at least one undefined bootstrap statistic, all for
precision: 1–713 undefined draws per affected contrast. They remain in the
family at p=1; no fraction cutoff rescues them. The other 64 contrasts have
numerically resolved decision gates at this fixed budget. Zero numerically
unresolved gates does not mean that the 16 untestable precision claims have
been validated. Paired sample counts include resampled multiplicities and vary
with eligibility; they are not counts of newly observed patients.

The useful development conclusion is narrower than the fixed-reference point
comparison: continuous amplitude estimation can improve correlation over some
weaker outputs, but this experiment does not establish a joint accuracy/error
advantage over the strongest retained baselines. Keep both estimators research
only. Do not retune a reference, comparator family, caller population or interval
rule to obtain a win. Further gains require a justified mechanism or additional
independent information, not a larger Monte Carlo budget for these already
resolved numerical decisions. The full manuscript/SOTA goal is not achieved.

All numerical outputs, including missingness and both signed/oriented point
effects, remain in the linked summary. The
[resource journal](../dev/bulk_amplitude_uncertainty_20261109/resource_summary_final.json),
[accounting](../dev/bulk_amplitude_uncertainty_20261109/accounting_final.tsv),
[audit inventory](../dev/bulk_amplitude_uncertainty_20261109/audit_inventory.json),
[specimen correspondence](../dev/bulk_amplitude_uncertainty_20261109/specimen_correspondence.json)
and [source integrity](../dev/bulk_amplitude_uncertainty_20261109/source_integrity_checkpoint.json)
are retained. Maximum process RSS across the batch is 1,883,396 KiB, below the
8GB allocation. The first and last draw chunks took 39.71 and 70.38 seconds;
accepted summary took 17.02 seconds and 898,704 KiB. Short-job scheduler RSS
can be zero; process high-water records supply the measurement.
All 75 production/native/test source files remain unchanged. The full production
suite was not rerun for this research-only analysis. No method is promoted and
no claim of biological uncertainty coverage or independent superiority is made.
