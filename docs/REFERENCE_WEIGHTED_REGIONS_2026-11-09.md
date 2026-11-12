# Reference and count weighting do not yet justify another bulk replacement

Reference-weighted native-region medians reach mean comparison Pearson .752184,
versus .738056 for the installed native-region option and .730261 for strong
pseudobulk. This is a promising point result, not established superiority. The
declared combined simultaneous-interval/Holm rule does not establish added
benefit over either retained native regions or the strongest baseline. Keep this
research candidate experimental. A matched count-only ablation reaches .739673;
its weighted mean falls to .411486 and is rejected. Preserve all outcomes and
comparators. No package, default, scientific threshold or frozen result changes.

## Measurement, estimand and competing explanations

The estimand is a regional location of programme-adjusted natural-log relative
RNA abundance against the specified normal panel. Median and mean are different
location functionals unless a common regional RNA effect is justified. Neither
is automatically DNA copy number. The fitted-bin PLN working model, gene
expectation allocation, half-count transformation, centring and native RNA-HMM
boundaries are inherited unchanged from the installed native-region option.

For each gene, estimate spread from the 20 leave-one-out normal-panel gene RNA
effects as MAD/normal-75th-quantile, squared. Add 1/(query fitted gene count+.5).
Use inverse square-root weights in the regional median or inverse weights in the
regional mean. This is a descriptive weighting heuristic, not an exact log-count
variance, confidence SE or calibrated sampling model. Reference measurement
noise may contribute twice. Between-normal spread may encode biological state,
reference mismatch or coverage differences as well as measurement reliability.

The count-only ablation sets the spread term to zero while preserving exactly the
same gene effects, boundaries, support, half-count and scoring. It distinguishes
count precision from learned reference variability as an explanation of the point
gain. The unweighted regional mean separates weighting from changing the location
functional. The reference-weighted median exceeds the count-only median by
.012511 in the point comparison, but the conditional simultaneous uncertainty
does not establish that mechanism. Reference variability must transport to query
specimens; robust aggregation alone cannot identify dosage against coordinated
transcriptional state, purity or reference-composition changes.

Corresponding independently processed DNA at matched resolution is the external
profile measurement. These existing specimens are development/regression evidence,
not unopened evaluation. Independent biological replication and appropriate truth
are still needed. No evidence here validates individual clone membership, spatial
assignment, absolute DNA dosage, purity, ploidy, interval coverage or the 5%
whole-specimen clone false-split target.

## Locked design, implementation and exact controls

The [reference design](../dev/reference_weighted_regions_20261109/DESIGN.md)
freezes two estimators and the matched unweighted mean before their runs. Source
counts, all 68 gene-effect/expectation contexts, the 20-normal panel and eight
previously fixed reference resamples come from the programme-gene experiment.
There are 24 distinct calibration and 24 distinct comparison tumour patients.
Only 23 comparison correlations and 13 loss recalls are finite. Normal candidate
rows are explicitly missing because normals supply the variance reference;
they are not silently included in tumour performance.

RNA predictions are saved before opening DNA scores. All 18 preceding method
scores reproduce exactly; adding three gives 21 methods. Exact dyadic integer
weight sums resolve the weighted-median half-mass condition without an epsilon.
Independent Fraction medians, adjacent-float half-mass examples, Fraction means,
gene permutation, weight scaling and chromosome/arm/gap resets pass. The original
native median also reproduces exactly. These are computational identities, not
scientific validation. Implementation:
[reference_weighted_regions.py](../benchmarks/hardening/reference_weighted_regions.py).

The [count-only design](../dev/count_weighted_regions_20261109/DESIGN.md) is written
after the reference uncertainty result and before count-only outcomes. It freezes
two count-only estimators, retains all 21 preceding methods exactly, and makes 23
methods total. Its additional Fraction fixture verifies removal of only the
reference-spread component. Both candidates are retained in uncertainty because
the median passes the declared positive-point-gain continuation gate against both
native median and strong pseudobulk. This continuation gate is not a promotion rule.

## All primary profile tradeoffs

| Method | Mean Pearson (23) | Loss recall (13) | Loss precision (calling patients) | Loss FPR (23) |
|---|---:|---:|---:|---:|
| Installed native-region median | .738056 | .495375 | .348000 (15) | .037392 |
| Reference-weighted median | .752184 | .577106 | .411623 (16) | .042304 |
| Reference-weighted mean | .748044 | .577106 | .407853 (16) | .047028 |
| Count-only weighted median | .739673 | .546658 | .384048 (17) | .043504 |
| Count-only weighted mean | .411486 | .072527 | .027192 (16) | .053892 |
| Unweighted regional mean | .669541 | .386401 | .241267 (14) | .066690 |
| Strong naive / own partition | .730261 | .399679 | .358571 (14) | .041244 |
| HD precision-weighted path mean | .733069 | .577106 | .458667 (15) | .046896 |

Reference median versus native median has 19 correlation wins/4 losses, and
versus strong naive 19/4. Count median versus native has 12/11 and versus strong
naive 13/10. Count mean loses all 23 correlation comparisons against native and
strong naive. Recall improvements coincide with increased FPR relative to native;
there is no claim that all metrics improve. Precision is conditional on emitting
calls, with unequal patient denominators, and is not an unconditional paired
success rate. Every patient, failure and comparator remains in each task's
bulk_v1/patient_scores.tsv, summary.tsv and report.json.

## Conditional uncertainty and the unchanged decision rule

The [reference uncertainty design](../dev/reference_weighted_regions_20261109/UNCERTAINTY_DESIGN.md)
freezes its analysis after primary point results and before resampling. Its family
has 24 contrasts: two candidates versus native median, old gene predictor, strong
naive and HD precision path mean, for correlation, loss recall and loss FPR.
The [count uncertainty design](../dev/count_weighted_regions_20261109/UNCERTAINTY_DESIGN.md)
uses another declared 24-contrast family, replacing old gene predictor with the
reference-weighted median. Both families remain recorded. Neither controls the
cumulative adaptive development history or supplies independent confirmation.

Each uses 20,000 nested draws: one of eight fixed normal-reference resamples,
resampled calibration patients and resampled comparison patients. Re-estimate the
same exact patient-balanced 5% neutral-window loss threshold in every draw. This
5% bulk window target is separate from the specimen clone error target. Reference
seeds are 2026110916/17; count-only seeds are 2026110918/19. Exact primary replay,
2,214 exact threshold checks and 5,184 direct score-lookup checks pass in EACH
analysis, with independent fsum checks and binomial Monte Carlo endpoint brackets.

Promotion requires BOTH the existing Holm-adjusted centred-bootstrap p < .05
and simultaneous basic interval lower endpoint > 0. Those summaries are not dual
under asymmetric bootstrap distributions. A favourable interval alone cannot
replace the combined frozen rule. Approximate conditional bootstrap behaviour
does not establish repeated-workflow coverage or calibrated biological errors.

| Pearson contrast | Paired gain | 24-family interval | Holm p | Added benefit established? |
|---|---:|---|---:|---|
| Reference median vs native median | +.014128 | [.004694,.031405] | .134193 | No |
| Reference median vs strong naive | +.021923 | [-.004493,.051521] | .451977 | No |
| Reference median vs HD precision | +.019116 | [-.003762,.038291] | .119594 | No |
| Reference median vs old gene predictor | +.071939 | [.021497,.124292] | .002400 | Yes, against an already surpassed comparator |
| Count median vs native median | +.001617 | [-.012689,.020477] | 1 | No |
| Count median vs reference median | -.012511 | [-.029869,.010606] | .743713 | No |
| Count median vs strong naive | +.009412 | [-.019547,.037725] | 1 | No |
| Count mean vs native median | -.326570 | [-.448148,-.094284] | .001200 | Harm supported conditionally |
| Count mean vs strong naive | -.318775 | [-.440727,-.085931] | .001200 | Harm supported conditionally |

Reference median's loss-recall gain versus native is .081731 with interval
[-.208193,.432317]; its FPR improvement is -.004912 with interval
[-.066118,.045776]. Count median's recall gain versus native is .051282 with
[-.268084,.255348]; FPR improvement is -.006112 with [-.073565,.028928]. These
loss claims remain unresolved. All 48 contrasts are saved, including failed and
harmful candidates. More resampling cannot create additional patients or erase
the lack of a verified unopened evaluation panel. CPTAC is already exposed in the
historical ledger and is not an unopened replacement; no new cohort is opened.

## Compute, attempts and preservation

All 24 jobs below complete 0:0 on CPU compute nodes using the same read-only
installed noneditable wheel, native binary and neutral-cwd import checks.

| Jobs | Role | Account/QOS; CPUs / memory / time |
|---|---|---|
|23149988|reference numerical checks|nslab; 1 / 4GB / 1h|
|23149989|reference primary|innovation; 4 / 16GB / 2h|
|23150140–23150148|reference fits 000–008|alternating innovation/nslab; 1 / 4GB / 30min|
|23150149|reference uncertainty|nslab; 1 / 8GB / 1h|
|23150221|count numerical checks|nslab; 1 / 4GB / 30min|
|23150222|count primary|innovation; 1 / 4GB / 30min|
|23150247–23150255|count fits 000–008|alternating innovation/nslab; 1 / 4GB / 30min|
|23150256|count uncertainty|nslab; 1 / 8GB / 1h|

Reference source_v1 runs numerics/primary; source_v2 is unexecuted; source_v3
runs refits/uncertainty after correcting design presentation. Count source_v1
runs numerics/primary; source_v2 adds the explicitly configured uncertainty family
and runs continuation. The reference-uncertainty defaults preserve its completed
configuration exactly. All snapshots, scripts, designs and transitive scorers
are hashed; exact command arrays are in each jobs.json. Resource/accounting and
source-integrity records are saved alongside outputs. No GPU or new dependency
is used. The previous installed-wheel 421 Python/4 Rust passes and one GPU-only
skip are separate evidence; that suite is not rerun for research-only code.

Wheel SHA256: 4e6e5f42ec26911efc9d118c511563c5cf1dfa2026a7456fe39d08765230dd07.
Native SHA256: e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7.
Shared installations and all frozen outputs remain unchanged. The new weighting
is not integrated because its incremental benefit gate has not passed. Keep the
installed native-region opt-in and reject the count-weighted mean.

[Reference results](../dev/reference_weighted_regions_20261109/bulk_v1/report.json),
[reference uncertainty](../dev/reference_weighted_regions_20261109/uncertainty_v1/report.json),
[count results](../dev/count_weighted_regions_20261109/bulk_v1/report.json),
[count uncertainty](../dev/count_weighted_regions_20261109/uncertainty_v1/report.json).
