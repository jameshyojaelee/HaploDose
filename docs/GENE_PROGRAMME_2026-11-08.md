# Direct gene aggregation with programme-effect transfer

Numerical qualification: the later [Monte Carlo precision audit and fixed-budget
refinement](PROFILE_MC_PRECISION_2026-11-08.md) first found unresolved tail precision,
then resolved it with one million draws and independent checks. The refined
programme-direct versus programme-median Pearson Holm p is0.0186340, with
numerical upper0.0234591. Biological/reference limitations remain; original
numbers below are unchanged.

The existing programme correction and direct gene aggregation do not demonstrate
an additive benefit in this matched experiment. Retain the simpler direct-gene
research reference; do not promote the transferred correction or change the
production bulk model from these results. The earlier programme ablation on the
summed-bin model remains separate evidence, not invalidated by this experiment.

## Fixed question and measurement

The experiment reuses nine frozen raw-STAR fits: 20 normal donors, 24 calibration
patients and 24 previously inspected comparison patients. It is development and
regression evidence, not independent evaluation. All 15 prior methods remain.

For each bin, the correction is the frozen programme-corrected centered log ratio
minus the corresponding uncorrected summed-count ratio. This same additive RNA
effect is transferred to every gene in that bin. Compare direct gene aggregation
and intermediate bin medians with this correction added. Reference, normalization,
centering, common support and the final RNA-only segmentation operator remain
fixed. Bin and gene half-count offsets differ; this is explicitly an effect-scale
diagnostic, not a gene-level count likelihood or reconstructed expected counts.

Nine exact replays recover the previous uncorrected direct-gene and median-bin
paths. Count conservation, sample alignment, bin correction algebra and within-bin
contrast preservation pass before scoring. Calibration-patient data alone set loss
thresholds; DNA does not select segmentation boundaries.

## Results

Six new contrasts across three endpoints join all 75 earlier endpoints in a
93-member family. The original seed and 20,000 nested patient/reference bootstrap
resamples are retained. Intervals below are simultaneous Bonferroni basic bootstrap
intervals; p-values use Holm93. They account for the finite eight-reference-refit
ensemble, whose adequacy and tail coverage are not independently established.

| Pearson contrast | Mean paired effect | Simultaneous interval | Holm p |
|---|---:|---:|---:|
| Transferred programme + direct genes minus uncorrected direct genes | −0.01569 | [−0.03719, 0.02459] | 1.0 |
| Transferred programme + direct genes minus transferred programme + intermediate bins | 0.03162 | [0.00737, 0.06252] | 0.01155 |
| Transferred programme + direct genes minus HaploDose under the same operator | 0.03765 | [−0.02324, 0.07692] | 1.0 |
| Transferred programme + direct genes minus naive under the same operator | −0.02117 | [−0.04372, 0.02276] | 1.0 |

The paired correlation population has 23 of 24 comparison patients; one remains
unevaluable. None of the new Spearman contrasts passes the family criterion.
Loss recall is evaluable in 13 patients; all 11 missing cases remain explicit.
Adding the correction to direct genes gives mean paired recall change 0.0 with
simultaneous interval [−0.3932, 0.3668]. Every new recall contrast is inconclusive.
An interval spanning zero is not evidence of equivalence or noninferiority.

Descriptively, the direct transferred path has mean loss recall 0.3997 and achieved
loss-side FPR 0.04337. Its uncorrected direct counterpart has the same mean recall
and FPR 0.04556. These FPR point estimates do not certify 5% error control. Mean
precision is 0.3124 versus 0.3569, each defined in 15 patients; this is descriptive,
not a demonstrated precision improvement. Classified-truth coverage stays 1.0
on the prespecified common support, not the entire genome.

## Decision and validation

[Observed pattern, not a guarantee] The intermediate-bin loss persists after this
matched programme transfer. The correction itself shows no task-specific benefit
on the direct-gene path, so the additional component is not retained as an improved
method. This does not show all gene-level programme models fail, nor justify
developing a more elaborate one without a new, testable mechanism. No superiority
over the simple baseline, count uncertainty calibration or biological CNA validity
is established.

Panel job 22972644 completes in 24m50s, scheduler peak RSS 643,140 KiB. Independent
audit 22972645 completes in 12s and verifies 99 input hashes, all 93 paired effects
and sample counts, the unchanged patient scores of all 15 prior comparators,
423 baseline files and 66 unchanged package/test/native files. All nine refit
arrays, operator details, scores, intervals and execution identities are retained
under `dev/gene_programme_20261108/gene_v1`; the audit is `audit_v1/report.json`.
No frozen artifacts, production defaults, comparator families or scientific
criteria were overwritten. The existing all-missing-window warnings retain missing
support rather than substituting a value.
