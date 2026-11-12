# Bulk programme attribution and error-bound design

The goal remains accurate regional dosage and appropriately supported CNA and
clone inference. This investigation tests whether removing reference-derived
expression programmes repairs the bulk profile gap. It does not: the component
helps some profiles, while the strongest matched segmentation comparison remains
unresolved. No production inference change or new validation claim is promoted.

## The tested mechanism

The estimand is a patient-average difference in DNA-profile association or
loss-detection performance between RNA profiles with and without programme mean
correction. The original 48 raw-STAR COAD/READ patients, 20 reference donors,
24 calibration/24 previously inspected reserved roles and eight donor-reference
bootstrap refits are unchanged. They remain development evidence, not a new
independent evaluation.

For each frozen fit, the new runner reconstructs raw bin counts using the
recorded gene map. Count conservation, sample/gene identity and recorded binned
totals must agree exactly. The no-programme ratio is
`log((Y + 0.5)/(E0 + 0.5))`, where `E0` is the sample's binned count total times
the pooled-normal bin rate; the result is centred by its bin median. This is the
RNA ratio produced by the existing `fit_sample` algebra with an empty programme
basis. No HMM state refit, new binning, new reference selection or DNA-informed
boundary fitting is performed.

The strongest competing explanations are removal of dosage-shaped signal versus
removal of transcriptional nuisance. The distinguishing observation is the paired
DNA agreement after matching genes, reference membership, callable support and
the smoothing operator. Raw ratios, the original 101-gene median, and the exact
RNA-only partition operator are evaluated separately. All units are natural-log
RNA effects. The existing HD segmented output remains a comparator; it is not
misrepresented as a newly fitted no-programme HMM result.

All seven previous methods and 24 comparisons are retained. The three
programme-on/off contrasts and the no-programme partition versus naive partition
contrast add 12 endpoints, giving a fixed 36-member family. Holm tests and
Bonferroni basic-bootstrap intervals use 10,000 nested reference/calibration/
reserved-patient resamples. Reference refits are not independent patients. Eight
reference draws and previously inspected patients cannot establish interval
coverage or prospective superiority.

## Results and decision

Effects below are programme-corrected minus no-programme profiles.

| Matched stage and endpoint | Mean paired effect | Simultaneous interval | Evaluable patients |
|---|---:|---:|---:|
| Raw ratio, Spearman | +0.0410 | [0.0154, 0.0675] | 23/24 |
| Median-smoothed ratio, Pearson | +0.0378 | [0.0128, 0.0596] | 23/24 |
| Own partition, Spearman | −0.0038 | [−0.0250, 0.0324] | 23/24 |
| Own partition, Pearson | −0.0028 | [−0.0270, 0.0633] | 23/24 |
| Own partition, loss recall | −0.0098 | [−0.3937, 0.2380] | 13/24 |

Both positive correlation results pass the complete multiplicity rule
(Holm p=0.00360). Own-partition differences do not establish a benefit or harm.
Removing programme correction therefore fails as a repair for the observed
profile gap; retain the existing component. These results do not show a general
advantage of the full package over all comparators.

The no-programme partition still has lower Pearson correlation than the naive
partition: difference −0.0560, simultaneous interval [−0.0921, −0.0125], Holm
p=0.00360. Gene aggregation and reference construction remain possible causes;
they were not intervened on here and are not established explanations.

Loss recall alone must not be read as a fixed-error advantage. Programme
correction raises the raw ratio's mean recall from 0.3089 to 0.5460, but achieved
mean loss-side FPR also changes from 0.0418 to 0.0505, and the simultaneous recall
interval includes zero. Median-smoothed FPRs are 0.0431 without programmes and
0.0429 with them. For own partitions they are 0.0356 and 0.0423. All methods retain
the same classified-truth coverage; missing endpoint denominators remain explicit.
There is no cutoff tuning, default removal or claim of universally improved loss
detection.

## Independent error-bound design calculation

A separate numerical job reads no simulation or biological outcomes. It audits
what the existing 5% upper-confidence-bound criterion can demonstrate with
200 evaluation specimens, under independent binomial sampling.

| Confidence family | Maximum acceptable errors / 200 | Upper bound at that count | Upper bound with one additional error |
|---|---:|---:|---:|
| One prespecified endpoint | 4 | 0.04518 | 0.05184 |
| 27 simultaneous endpoints | 1 | 0.04194 | 0.05134 |
| Current 45 simultaneous endpoints | 1 | 0.04467 | 0.05427 |

At a true 5% error rate, the probability of meeting the current simultaneous
criterion is 0.000404, or 0.0404%. At true rates of 1% and 0.5%, it is 0.4046 and
0.7358 respectively. A hypothetical procedure would need a true rate at most
approximately 0.00412 for an 80% probability of passing this fixed 200-specimen
criterion. These are prospective operating-rate scenarios, not estimates of
HaploDose's performance.

Independent binomial-tail inversion agrees with beta-quantile bounds to within
2.64e-16. The calculation also verifies the original requirement of at least
59 independent null specimens for zero errors to satisfy the single-endpoint
one-sided 95% bound. Simulated seeds cannot supply biological replication.

Failing to demonstrate an upper bound at or below 5% does not itself establish
that the true rate is above 5%. That latter conclusion would require contrary
evidence, such as a lower bound exceeding the target. The pass criterion remains
unchanged. No extra evaluation replicates or new operating cutoff is selected
to rescue the running experiment.

The report includes a mathematical calibration order-statistic calculation as
a future design reference. Its guarantees require exchangeable independent null
scores; they do not protect against biological misspecification or certify
membership. It is not implemented as an inference rule or applied to current
outcomes. Lowering an operating error rate can reduce power, which must remain
part of any later test.

## Verification and evidence

`dev/programme_attribution_20261108` contains the locked design, dirty baseline,
executed source snapshots, source/input hashes and new outputs. Job22962704
completed the component comparison in 8m59s; audit22962740 verified 106 input
hashes, 378 baseline files, 65 unchanged package/test/native files, exact
preservation of the seven original comparator results and all 36 paired effects.
Job22962859 completed the independent error-bound calculation. All executed on
CPU compute nodes with the existing isolated validated wheel.

The new scripts are `benchmarks/hardening/programme_operator_ablation.py`,
`programme_operator_audit.py` and `error_control_design.py`. Production inference,
native code, frozen results and original tests are unchanged. The full package
suite was not rerun because no production code changed. README now states the
existing bulk functionality and links the current evidence report.

Computational checks support the transformations and paired calculations.
Identification remains conditional on reference composition, measurement-scale
assumptions and the specified DNA/RNA correspondence. Biological validity,
routine interval coverage and competitor superiority remain unestablished.
The earlier exact-confirmation false-split/power jobs are still running at this
checkpoint; no partial outcomes have been used to change their method or criteria.
