# A fixed median reference ensemble does not restore bulk loss recovery

The implemented median ensemble fails the frozen performance gate. Mean
Pearson is0.770177 versus the retained arithmetic ensemble0.770411; loss recall
stays0.540476 rather than recovering the preceding single-reference0.602747.
Reject this correction and stop its uncertainty expansion. All84 old method
scores reproduce EXACTLY; both compute jobs complete0:0. No default changes.

## Intervention and the result that would have justified it

The preceding arithmetic reference ensemble improves profile accuracy but
attenuates loss recovery. Test whether a minority of reference fits washes out
otherwise recurring RNA effects:take the fixed median of ALL8 existing
reference fits001..008, with the midpoint of middle order statistics. Apply
the IDENTICAL operator to all39 eligible methods, including strong pseudobulk.
Retain ALL84 preceding methods for123total. No fitted weight, trimming fraction,
reference selection or extra resample. This is one prespecified intervention,
not a choice among medians selected after DNA scoring.

Estimand:regional relative RNA abundance at the existing validation-window
resolution under those empirical reference perturbations. The working model
is robust aggregation of the retained RNA predictions. The identifying
restriction is that most reference fits preserve the relevant signal. Shared
lineage mismatch or regional transcription can survive all8 fits and resemble
DNA dosage. Independent corresponding DNA at matched resolution, plus new RNA
references, would distinguish genomic prediction from shared reference bias.
No absolute DNA state, clone membership or spatial assignment follows.

The smallest decisive comparison uses the same48 patients,24/24 split,
23 informative evaluation correlations,13 loss-recall patients, features,
callable support, truth, scorer, threshold-calibration rule and failure handling.
The RNA-only job locks all39 outputs before a dependent scorer accesses DNA.
Independent Fraction midpoint, permutation and minority-outlier fixtures pass.
All eight reference supports and specimen orders must match exactly.
[Design frozen before outcomes](../dev/robust_reference_profile_20261110/DESIGN.md).

Primary must beat ALL122 alternatives in paired mean Pearson. The joint point
gate additionally requires no mean paired loss-recall decrease or FPR increase
against retained single-reference `averaged_gaussian_reference_median`.
These are direction-of-benefit gates, not changed alpha/error thresholds.
A point pass would only permit a separately frozen uncertainty continuation;
the earlier failed uncertainty and independent-validation claims remain.

## Complete performance and decision

| Method | Mean Pearson /23 | Loss recall /13 | Neutral loss FPR /23 |
|---|---:|---:|---:|
|New fixed median ensemble|0.7701774107|0.5404762|0.0412094|
|Retained arithmetic ensemble|0.7704112303|0.5404762|0.0447602|
|Retained single-reference boundary average|0.7648167082|0.6027473|0.0437372|

Primary-minus-arithmetic Pearson is-0.0002338197,11 wins/12 losses. This is
the only negative primary Pearson contrast in the full122 family, but it
decisively fails the unchanged all-comparator point rule. Primary-minus-
single-reference paired recall is-0.0622710623; FPR is-0.0025277663. The FPR
point reduction does not justify a loss-recovery claim or override the joint
gate. No new uncertainty study, selected secondary or default promotion.
[All123 methods and122 contrasts](../dev/robust_reference_profile_20261110/score_v1/report.json),
[every patient score](../dev/robust_reference_profile_20261110/score_v1/patient_scores.tsv).

Computational assessment:all required identity/support fixtures and84old-score
exact controls pass. Statistical assessment:fixed conditional regression point
gates fail; no new coverage/superiority result. Biological assessment:no new
replication, corresponding DNA or membership truth; the exposed cohort remains
regression-only. The result fails to support minority-reference outliers as a
sufficient correction for the observed recall loss. It does not prove that
all reference effects or all other mechanisms are absent.

## Exact code, execution and preservation

[Implementation](../benchmarks/hardening/robust_reference_profile.py),
[locked RNA arrays/hash](../dev/robust_reference_profile_20261110/aggregate_v1/RNA_LOCKED.json).
Both jobs use CPU/v6&amd on005, capturedruntime003, the already validated
noneditable wheel and an immutable source snapshot. No package, native,
dependency, original test, criterion, old source or frozen result is changed.
No full package/GPU suite rerun was needed for this research-only experiment.

| Job | Role | Account/QOS | CPU / RAM / requested time | Outcome |
|---|---|---|---|---|
|23158600|39 RNA medians|innovation|1 /4GB /1h|complete0:0,12s|
|23158601|123-method dependent scoring|nslab|1 /4GB /1h|complete0:0,6s|

[Commands](../dev/robust_reference_profile_20261110/jobs.json),
[accounting](../dev/robust_reference_profile_20261110/accounting_final.tsv),
[resources](../dev/robust_reference_profile_20261110/resource_summary_final.json),
[preservation](../dev/robust_reference_profile_20261110/preservation_final.json).
The clone recovery/full-workflow false-split work continues independently;
this negative result does not alter its gates or establish the full goal.
