# Membership weighting under coherent gene splits: partial recovery

Weighting both confirmation and merge recovers4/5 original cases, versus3/5
for the same fixed known-gene-noise splits with neither or either stage alone.
20261110 is recovered;20261109 still fails. Production remains2/5. This is a
research interaction, not a sufficient correction, a production change or a
biological validation claim.

## Estimand, intervention and competing explanation

The estimand is recoverable regional RNA-profile populations under the existing
working NB/pooled-profile models. The previous observational trace localized
20261110 loss before matched confirmation and20261109 loss at merge. The
locked `dev/gene_thinning_weighted_20261109/DESIGN.md` tests the pre-existing
membership_weighted_stages implementation, without a new fitted constant.
Confirmation weights hard proposal groups by fitted membership margins, uses
squared weights for conditional sampling variance, then performs matched support.
Merge weights use existing proposal membership margins. The same gene test half
is supplied to the veto; full-count profiles/nuisance and stopping remain.

These are changed estimators, not a claim that historical hard pooling was a
coding error. The confirmation arm replaces initial hard-pool agglomeration with
the existing weighted comparison. Numeric recovery assertions and existing
posterior/bin/gene/cell support requirements are unchanged. Fitted weights are
selected: conditional variance formulas omit their estimation and selection
uncertainty. Known source gene size10 remains unavailable oracle information.
Continuous expression state, reference error and full-count nuisance leakage
remain alternatives to a DNA-defined distinction.

| Seed | Gene hard accuracy | Confirmation only | Merge only | Both | Both passes all original criteria |
|---|---:|---:|---:|---:|---|
|20261031|0.836032|0.836032|0.836032|0.836032|yes|
|20261107|0.738095|0.738095|0.738095|0.738095|yes|
|20261108|0.790514|0.790514|0.790514|0.790514|yes|
|20261109|0.601202|0.601202|0.601202|0.601202|no|
|20261110|0.599602|0.598000|0.599602|0.810000|yes|

All four arms' accuracy uses the original tumour-call-conditional criterion and
requires separate exactly-two-group, normal-FPR, tumour-TPR, one-group-S2 and
chr12-only support checks. In the combined20261110 arm, normal false calls6/90,
tumour detection500/510, exactly two S1 tumour groups, one S2 group and only chr12
distinctions pass.20261109 retains one tumour group, accuracy0.601202 and no
chr12 distinction; its failures remain. No significance or population performance
estimate follows from five previously exposed mechanism cases.

## Independent checks and decision

Independent audit_v2 recounts all90 original criterion values including failures,
with explicit hashes for fitted inputs, original truth, stage traces and split
arrays. Every array in all45 saved fold artifacts equals the corresponding
original gene arm, including actual veto subsets and count conservation. This
removes changed random splits as an explanation of this within-gene stage
comparison; it does not resolve the earlier bin-versus-gene comparison's changed
draw realizations. Actual stage traces verify both weighting mechanisms execute.
The same scorer compares every original criterion; no failed method disappears.

The three inference jobs23119669/670/671 complete0:0, CPU4/8GB/2h requests on
innovation/nslab. Audit23119726 failed because its stage-trace path omitted the
panel directory. Only that path is corrected; failed audit_v1 and source stay
frozen. Audit23119761 completes0:0 with the same criteria and fold checks,
CPU1/8GB/1h nslab. Five attempts: four successful, one preserved audit failure.
No inference rerun. Peak process RSS1616784KiB, below1.55GiB; actual timing,
memory, backend, wheel/native identities and commands are in the task records.

Research changes add an optional stage selector to gene_thinning_recovery.py
and the independent gene_thinning_weighted_audit.py. Executed sources are frozen
in source_v1/audit_source_v1/v2. Per-arm manifests retain the foundational
gene-thinning design; the new interaction design and all arm inputs are bound by
the independent audit manifest and recorded submission commands. All76 protected
package/native/test/README/build files match the initial dirty baseline. The
installed wheel's395 Python/4 Rust passes and1GPU-only skip are previous checks;
full tests were not rerun for this research-only intervention.

Reject this as a sufficient recovery correction or default. The locked design
does not warrant a larger calibration sweep for arms that fail original recovery.
No new complete-workflow null/power, repeated-workflow coverage or biological
test was run. The5% whole-specimen upper confidence bound and simultaneous
detection-power requirements remain unchanged and unachieved; biological n=0.
A calibrated biological result would require independent nuisance/confirmation
evidence, repeated-workflow error and coverage tests, and corresponding DNA with
documented sensitivity. More repetitions of these five simulations cannot
substitute for that missing evidence.
