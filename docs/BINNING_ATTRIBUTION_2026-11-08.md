# Intermediate binning contributes to the bulk profile gap

The objective is to distinguish added genomic information from differences in
RNA preprocessing and smoothing. This development comparison finds a loss of
Pearson agreement when the naive profile undergoes the same intermediate gene
binning used in the preceding HaploDose median-aggregation diagnostic. It does
not establish that all aggregation is harmful, that the residual methods are
equivalent, or that a revised package has passed independent evaluation.

## Locked intervention and estimand

The added method starts from the exact frozen naive gene effects: full-library
CPM+1, the mean normal-donor log abundance, and the original centering universe.
It takes the median within each frozen HaploDose bin, expands those values back
to their constituent genes, then applies the same evaluation-window median and
RNA-selected partition operator. All mapped RNA genes contribute to bins before
the common DNA-mapping mask is applied. This changes only the intermediate
aggregation. Naive gene effects reproduce exactly in all nine frozen fits.

The estimand is the biological-unit mean paired difference in RNA–DNA agreement
or loss recall on the declared common support. The working RNA effect is a
relative log-abundance contrast, not absolute DNA copy number. An alternative to
a better dosage model is altered noise suppression or signal loss from nested
aggregation; matching that operation distinguishes these explanations here.
Medians do not commute, as verified by an independently calculable six-gene
fixture. Identity permutation and translation invariance also pass.

The same48 patients retain24 calibration and24 previously inspected comparison
roles. These are exposed development/regression data, not a new evaluation.
Twenty normal donors and eight frozen reference resamples are reused. All13
earlier methods and54 endpoints remain; one new method and nine endpoint
contrasts produce a63-member family. The procedure retains20,000 nested
reference/calibration-patient/comparison-patient bootstrap draws, Holm tests and
Bonferroni basic intervals. Eight reference draws and extreme bootstrap tails do
not establish nominal interval coverage. DNA segment uncertainty is not resolved
by this comparison.

## Results and decision

| Contrast | Endpoint | Mean paired difference | Simultaneous interval | Holm p |
|---|---|---:|---|---:|
| Binned naive − original naive partition | Pearson | −0.03128 | [−0.05057, −0.01712] | 0.00315 |
| Binned naive − original naive partition | Spearman | −0.01669 | [−0.03216, 0.00675] | 0.258 |
| Pooled-reference median − binned naive | Pearson | −0.01012 | [−0.01994, 0.00440] | 0.498 |
| HaploDose ratio partition − binned naive | Pearson | −0.02753 | [−0.06714, 0.05438] | 1 |

Correlation contrasts have23 evaluable patients; the24th remains missing rather
than disappearing from the population. Loss-recall contrasts have13 evaluable
patients and11 without an estimable loss-recall denominator. The binned-naive
recall difference is−0.10256, interval[−0.34799,0.29487]; it is inconclusive.
All nine newly added effects, missing counts and intervals remain in the output.

[Inference] Intermediate aggregation contributes to the Pearson gap on this
panel. The previously significant pooled-median versus naive gap becomes
inconclusive after matching binning; this is not an equivalence result. The
smallest justified follow-up is to bypass intermediate bins in the existing
pooled-reference gene-effect calculation, keeping its normalization, half-count
offset, common support and final RNA-only partition fixed. An improvement that
vanishes under those controls would reject bin bypass as a useful correction.
No production default is promoted from this result.

## Executed evidence

Job22964738 completed on an Intel v4 CPU node in22m52s, with4 CPUs and16GB
requested; process peak RSS was681,652KiB. Audit22964768 completed in11s with
2 CPUs and4GB requested, process peak RSS307,552KiB. The audit verifies89 input
hashes, all63 paired effects and counts, all13 retained comparators,398 baseline
files and65 unchanged package/test/native files. Original outcomes and source
snapshots remain untouched. Audit scope does not include bootstrap-coverage proof.

Evidence is under `dev/binning_attribution_20261108`, including
`binning_v1/comparisons.tsv`, `binning_v1/patient_scores.tsv` and
`audit_v1/report.json`. Analysis and scorer sources are bound by immutable
execution manifests; no production dependencies or scientific thresholds change.
