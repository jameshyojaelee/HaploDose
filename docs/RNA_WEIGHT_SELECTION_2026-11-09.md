# Held-gene RNA risk does not retain the reference-weighted profile gain

The implemented RNA-only selector chooses uniform native-region medians for
47/48 tumour samples and reference weighting for one. Every comparison-patient
profile is the unchanged uniform profile: mean Pearson .738056 versus .752184
always-reference-weighted and .730261 strong pseudobulk. Reject this selector at
its prespecified incremental-performance gate. Its larger reference-resampling
and uncertainty study is not run. This is a negative performance experiment,
not another package feature or a passing superiority claim.

## Proposed correction and scientific contract

Reference weighting gives a promising correlation point result but fails its
existing combined interval/Holm promotion criteria. Rather than select weighting
from DNA outcomes, the [frozen design](../dev/rna_weight_selection_20261109/DESIGN.md)
tests one global choice per sample using held-gene absolute RNA prediction error.
No new segmentation, reference, variance estimate or fitted DNA parameter enters.
The two candidate profiles are already fixed uniform and reference-weighted
medians on identical native regions. Every final selected profile must reproduce
exactly one of those two retained profiles.

The estimand is regional relative RNA location against the specified normal
panel. Ten folds, matching the existing RNA selector's fixed fold count, are
assigned by literal gene-ID rank modulo10 within each native region. Each held-out
gene is predicted from other genes in its region using each median. Sum unweighted
absolute errors over all eligible genes with math.fsum in fixed identity order.
Singleton regions lack a training gene and contribute to neither risk; all genes
remain in the final descriptive profile. The smaller risk wins globally for the
sample; overlapping numerical risk orderings select the simpler uniform profile.
The error budget is derived from operation count, binary64 precision and observed
magnitude, not a fitted scientific constant.

The identifying working assumption is that predicting unweighted regional RNA
variation is a useful surrogate for predicting matched-resolution DNA dosage.
Gene-specific transcriptional response, reference composition and expression
state are the strongest alternatives. A method can predict RNA better while
retaining transcriptional signal that does not report DNA. The corresponding
independently processed DNA supplies the external check of this surrogate.
The current exposed specimens are regression evidence, not unopened evaluation;
biological replication is unchanged.

This is explicitly conditional held-gene prediction: features, native regions,
gene expectations, programme correction and reference weights were fitted from
full RNA observations. The held-out effects do not enter the immediate training
median, but this is not independent confirmation or a claim that all nuisance
fitting excludes those counts. No calibrated confidence, genomic clone or absolute
copy-number inference follows from a low prediction error.

## Computational gates and fair comparison

[Implementation](../benchmarks/hardening/rna_weight_selection.py) passes prior
Fraction weight/median/mean fixtures and new analytically solvable RNA-risk
fixtures. In a ten-gene binary example, leave-one-gene uniform medians have exact
total error100; fixed weighted medians have error50. Independent scalar medians
agree. Altering only a held-out fold leaves every training prediction identical.
Gene permutation preserves folds by identity, both risks, final choice and final
profile exactly. Constant effects select uniform; gap resets isolate regions;
duplicate identities fail explicitly. The checks pass on compute nodes.

All48 tumour predictions, risk sums, selected estimators and numerical budgets
are written before DNA scoring in RNA_DECISIONS_LOCKED.json. Normal candidate
rows remain explicitly missing. Every selected profile equals its frozen full-data
source profile exactly. All23 preceding methods are retained, giving24 total;
their patient scores reproduce exactly using the same scorer, genes, samples,
callable windows and calibration threshold procedure.

There are24 distinct calibration and24 comparison tumour patients, with23 finite
comparison correlations and13 loss recalls. Precision is conditional on emitting
calls and has unequal denominators. Windows and cells are not biological
replicates. The patient-balanced5% neutral-window loss calibration is unchanged
and is distinct from the5% whole-specimen clone false-split target.

## All main performance tradeoffs and decision

| Method | Pearson (23) | Loss recall (13) | Loss precision (calling patients) | Loss FPR (23) |
|---|---:|---:|---:|---:|
| RNA-only selected median | .738056 | .495375 | .348000 (15) | .037392 |
| Uniform native median | .738056 | .495375 | .348000 (15) | .037392 |
| Always-reference-weighted median | .752184 | .577106 | .411623 (16) | .042304 |
| Strong naive / own partition | .730261 | .399679 | .358571 (14) | .041244 |
| HD precision path mean | .733069 | .577106 | .458667 (15) | .046896 |

Compared with the always-weighted median, the selector has four correlation wins
and19 losses, mean change -.014128. It has23 exact correlation ties against uniform
native median. The numerical support criterion is not relaxed to keep the
weighted option. This shows that this particular unweighted held-gene RNA-risk
rule does not retain the DNA correlation point advantage of reference weighting.
It does not establish that every RNA-only selection method fails, or that reference
weighting is itself validated for DNA. Its earlier uncertainty failures remain.

The predeclared continuation requires positive paired correlation gain over
uniform native median, always-weighted median AND strong naive, as well as using
both estimators somewhere. The latter is met over the48 tumours, but neither
incremental correlation requirement is met. Reject the addition and stop. The
planned eight reference resamples, exact primary replay and12-family nested
bootstrap are **not run**. Do not substitute the comparison-patient native profile
for the candidate that failed, change the loss after seeing DNA, or claim added
benefit from beating only the weaker pseudobulk point estimate.

## Execution, preservation and remaining evidence

| Job | Role | Allocation | Elapsed / exit |
|---|---|---|---|
|23150441|independent numerical fixtures|nslab CPU1/4GB/30min|18s,0:0|
|23150442|dependent primary performance|innovation CPU1/4GB/30min|12s,0:0|

Immutable source_v1 and initial dirty baseline are preserved. Both jobs use the
validated read-only noneditable native-region wheel, neutral-cwd imports and
recorded wheel/native/source identity. All actual process memory/runtime,
backend, environment, dependency and source-integrity records are saved under
dev/rna_weight_selection_20261109. Full package tests/GPU parity are not rerun for
this research-only addition; the earlier421 Python/4 Rust passes and one GPU-only
skip remain separate evidence. No production/native/original test, dependency,
shared installation, scientific criterion or frozen result changes.

Computational status: numerical identities, exact source-profile selection and
all23 previous scores pass. Statistical status: the proposed surrogate does not
deliver the declared incremental DNA performance; no new uncertainty claim.
Biological status: no new matched DNA, independent replication or unopened
validation panel. Clone recovery, specimen false-split/power and calibrated
uncertainty requirements remain unresolved. Further biological guarantees require
suitable independent measurements, not more draws from these same patients.

[Primary results](../dev/rna_weight_selection_20261109/bulk_v1/report.json),
[all patient scores](../dev/rna_weight_selection_20261109/bulk_v1/patient_scores.tsv),
[RNA decisions](../dev/rna_weight_selection_20261109/bulk_v1/RNA_DECISIONS_LOCKED.json).
