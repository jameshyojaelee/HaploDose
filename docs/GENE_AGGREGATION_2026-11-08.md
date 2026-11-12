# Direct gene effects recover information lost by intermediate medians

The goal is to improve bulk relative RNA profiles by retaining measured signal,
while separating preprocessing effects from additional genomic information.
Bypassing intermediate bins improves Pearson agreement within the existing
pooled-reference gene-effect calculation on the exposed development panel.
It does not establish superiority over the simple naive baseline or validate
count-model uncertainty, DNA causation or an independent evaluation claim.

## Intervention and identifying assumptions

The preceding [binning control](BINNING_ATTRIBUTION_2026-11-08.md) showed a loss
when intermediate bins were added to the naive baseline. This follow-up removes
that same operation from the pooled-reference median path: gene log effects go
directly to the fixed evaluation-window median and RNA-only partition. Mapped-gene
library depth, pooled normal gene rates, gene half-count offset and the original
bin-median centering scalar remain unchanged. The entire previous pooled-median
partition reproduces exactly in all nine fits before this operation is bypassed.

The estimand is the patient-level mean paired difference in RNA–DNA agreement or
loss recall on common callable support. Relative gene effects remain normalized
RNA measurements. The competing explanation for improved agreement is reduced
signal loss under nested medians, rather than a better reference or an added CNA
model. Holding those quantities and the final segmentation operator fixed isolates
that mechanism here. Median noncommutation, translation, within-bin permutation
and an independently calculated scalar gene-effect fixture pass.

This is a programme-free attribution path. It does not authorize removing
programme correction from HaploDose: the separate programme ablation found a
benefit. It also does not identify absolute copy number, purity or ploidy.

## Complete paired results

The same48 patients retain24 calibration and24 previously inspected comparison
roles, with20 normal donors and eight frozen reference resamples. There is no new
evaluation cohort. All14 previous methods and63 endpoints are retained; four
new contrasts across three endpoints yield75 simultaneous comparisons.
The procedure uses20,000 nested bootstrap draws, Holm adjustment and Bonferroni
basic intervals. Biological-unit uncertainty, missing denominators and common
support remain explicit; bootstrap coverage and DNA-truth uncertainty are not
certified by this analysis.

| Direct gene path minus comparator | Endpoint | Paired effect | Simultaneous interval | Holm p |
|---|---|---:|---|---:|
| Pooled-reference median bins | Pearson | 0.03592 | [0.01670, 0.05666] | 0.00375 |
| Binned naive | Pearson | 0.02581 | [0.01253, 0.04483] | 0.00375 |
| Original naive partition | Pearson | −0.00548 | [−0.01536, 0.00366] | 1 |
| HaploDose ratio partition | Pearson | 0.05334 | [−0.02261, 0.09698] | 1 |
| Pooled-reference median bins | Spearman | 0.02064 | [−0.00175, 0.03904] | 0.0667 |
| Pooled-reference median bins | Loss recall | 0.08571 | [−0.27159, 0.37510] | 1 |

Correlation contrasts have23 evaluable patients and one retained missing outcome.
Loss recall has13 evaluable patients and11 without an estimable denominator.
The direct path and original naive partition have equal mean recall in this
sample; their paired interval[−0.30278,0.25] does not establish equivalence.
All12 new endpoint contrasts are retained in the comparison table.

[Inference] The controlled Pearson improvement supports intermediate median
aggregation as a source of information loss in this panel. The original naive
baseline remains competitive. Retain the direct gene path as a simpler research
reference; no production default, calibrated confidence claim or superiority
claim is promoted. Any package integration must retain the established programme
controls and be tested at its actual exported resolution. Independent biological
validation remains necessary after a candidate is frozen.

## Reproducible checks

Job22967913 completes in24m48s on an Intel v4 CPU node with4 CPUs and8GB requested.
Independent audit22967914 completes in10s with2 CPUs and4GB. It verifies90 input
hashes, all75 paired effects and sample counts, all14 old comparators,408 preserved
baseline files and65 unchanged package/test/native files. Historical results,
original thresholds and selected specimen identities are unchanged.

Evidence is under `dev/gene_aggregation_20261108`, including
`gene_v1/comparisons.tsv`, `gene_v1/patient_scores.tsv`, all nine control
reconstructions and `audit_v1/report.json`. Candidate v1 was not executed;
v2 retained float64 gene values until the original window export cast, avoiding
an unintended additional precision change. No new production dependency was added.
