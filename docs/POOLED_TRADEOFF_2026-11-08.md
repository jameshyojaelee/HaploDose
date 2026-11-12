# Pooling initialization does not solve clone error and power

## Completed precision-corrected replay

The separate corrected-wheel998-fit replay and audit22972815 now complete. All
20 inference jobs use precision wheel
`60fe42a9eabcf9f30a154ea512d5723ee11136b88d230030be136a73f2d7fc15`.
This version predates the terminal stopping correction and does not validate the
current terminal-corrected package or the later joint-membership candidate.

The original state-calibrated gate retains15/200 false distinctions, with
simultaneous upper0.155379; detection is142/200, with lower0.595457. At0.005,
false distinctions are0/200 with upper0.037243, but detection remains18/200 with
lower0.036953, and exactly two groups are accepted in1/200. Thresholds0.01/0.025
now have6/200 state-null distinctions, upper0.092143, and respectively92/108
detected alternatives. All bounds belong to the full99-endpoint family.

Relative to the old wheel, the original gate gains nine alternative detections;
the27-comparison Holm p-value is0.105469, so no paired improvement is established.
Raw state-null distinctions remain85/200 and raw alternative detections199/200.
The same strong error/power/membership limitation persists; no default promotion.

Audit22972815 verifies998 paired inputs and their assignment files,5400 decisions,
99 independently inverted binomial bounds/pass flags,27 paired discordances,
nine independently matched membership endpoints, all20 wheel identities and420
frozen baseline files. Evidence is in
`dev/precision_null_20261108/{summary_v1,audit_v1}`. Old72-endpoint results below
remain historical results with their own source and confidence-family identities.

The complete old-wheel experiment and independent audit have finished. Pooling
initial tumour groups reduces raw correlated-state subdivisions relative to the
earlier path, but does not meet the unchanged whole-specimen error criterion at
the original operating point. A stricter gate meets that conditional simulation
bound while losing most detection and almost all exact two-group recovery.
No inference default or biological validation claim is promoted.

## Fixed experiment and estimand

The experiment uses998 independent simulation draws:199 calibration draws for
each of the matched-NB and correlated-state nulls, plus200 evaluation draws in
each null and the two-CNA-population condition. Every specimen has600 cells and
240 bins, with20-bin chromosome blocks. The complete binned initialization,
proposal search, fitting, merging, splitting and stopping procedure is rerun.
Reference and nuisance laws are fixed; public gene preprocessing, gene/allele
channels and biological replication are outside this experiment.

The primary event is any accepted additional final RNA-group distinction in a
one-CNA population. Detection means any accepted multigroup distinction in the
two-CNA condition; it does not certify correct group count or membership. All
raw groups, abstentions and one-to-one membership errors remain reported.
Continuous regional expression state is the principal competing explanation.

The original rank threshold0.05 remains primary. Thresholds0.005,0.01 and0.025
were fixed before these completed outcomes and assessed together, not substituted
for a failed primary result. All historical endpoints are retained in72
simultaneous one-sided bounds. Bounds condition on the realized calibration
bank; this exposed development experiment is not independent biological
evaluation or a universal calibration guarantee.

| State-calibrated rank threshold | False distinctions /200 state-null specimens | Simultaneous upper bound | Detection /200 alternatives | Simultaneous detection lower bound |
|---|---:|---:|---:|---:|
| 0.05, original | 15 | 0.15286 | 133 | 0.55130 |
| 0.025 | 5 | 0.08230 | 108 | 0.42503 |
| 0.01 | 5 | 0.08230 | 92 | 0.34793 |
| 0.005 | 0 | 0.03571 | 18 | 0.03799 |

The unchanged error requirement is upper bound<=0.05. Only the strictest row
meets it for the correlated-state condition. At that point, exactly two groups
are accepted in1/200 alternatives; mean retained-cell coverage is0.08793 and
mean correct one-to-one cell fraction over all cells is0.04098. Rejecting most
distinctions is not successful clone recovery. At the original threshold, the
133 detected alternatives likewise do not imply correct memberships: only27/200
raw alternatives have exactly two eligible groups, with mean raw ARI0.08119.

The NB-calibrated gate accepts85/200 state-null distinctions even at0.005,
demonstrating observation-model mismatch that stricter ranks within that bank
cannot resolve. All matched-NB null rates are0/200. Good behavior under the fitted
null does not establish validity under correlated cell state.

## Execution identity and checks

Inference jobs22964075/077/078/079/080 use the preserved OLD wheel
`dev/confirmation_20261107/environment`, SHA256
`93205a319db1c1000148bf82efa8dd4f14a7c60f492d73e1a48c485123f833e9`.
These results must not be attributed to the later graph-precision correction.
Summary22964088, audit22969949 and tradeoff22969969 complete successfully.
The audit verifies998 paired generated inputs and assignment outputs,1800
decisions, independent exact-binomial bounds, execution identities and bounded
source preservation. The tradeoff independently solves binomial-tail equations
and retains all72 endpoints. Historical outputs are unchanged.

Evidence is in `dev/pooled_initial_null_20261108/{summary_v1,audit_v1}` and
`dev/pooled_tradeoff_20261108/tradeoff_v1`. A separate corrected-wheel replay is
still active; its99-endpoint family will be reported independently. A newly
identified terminal merge-check issue is also being investigated. Neither a
retrospective error reduction nor extra computation closes the missing suitable
DNA truth and biological-null replication requirements.
