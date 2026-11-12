# Frozen same-source test of the learned DNA-response profile

The scientific goal is to determine whether the measured-DNA prediction gain
survives a frozen test on the 50 previously reserved patient-labelled units.
The full HaploDose objective remains genomic profile accuracy and supported
clone/spatial inference across modalities; this bounded test addresses one
profile mechanism and cannot establish that full objective.

## Completed development-null experiment

The [preceding development report](DNA_RNA_RESPONSE_DEVELOPMENT_2026-11-11.md)
retains the failed unbounded predictor, its traced RNA extrapolation mechanism,
the corrected training-support rule, every baseline and numerical check. The
corrected gene-only model is simpler and slightly better than the blend. The
blended model's extra fitting step has no demonstrated benefit.

All **99** fixed conditional permutations and their independent readers complete.
Writers 23225809, 23225810 and 23225811 finish 0:0 in 8m52s, 8m47s and 8m57s.
Summary 23225812 finishes 0:0 in 13 seconds. Each draw repeats reference fitting,
RNA support estimation, response slopes, inner weight selection, prediction
locking and scoring. Each reader checks 96 SciPy coefficient references and
independently reconstructs regional losses: **9,504** coefficient checks across
the null draws, in addition to the observed studies.

RNA rows are permuted within lineage, RNA-library and DNA-assay strata. There
are **97 movable** and **22 singleton** development units; singleton pairings
remain fixed. The permutation space has
46,167,057,931,450,682,815,272,714,240,000,000 elements. These draws reuse the
same 119 patient-labelled development units; they are not 99 new biological
specimens. Conditional exchangeability remains an identifying assumption.

| Fixed target | Observed minimum MSE gain over three baselines | Conservative maxT exceedances | Randomization p | Simultaneous 95% Monte Carlo tail upper |
|---|---:|---:|---:|---:|
| Original blend | 0.01111281288 | 0/99 | 0.01 | 0.03657574498 |
| Simpler gene-only | 0.01174918364 | 0/99 | 0.01 | 0.03657574498 |

The two-target max statistic and Bonferroni-protected exact Monte Carlo bounds
pass the declared conditional continuation rule. This is evidence of RNA/DNA
association for the fixed development workflow. It does not replay historical
method selection, prove independent biological prediction, validate DNA truth,
or establish the 5% whole-specimen clone false-split target. The clone error/power
and uncertainty failures remain unchanged.

[Frozen null protocol](../dev/dosage_response_null_20261111/DESIGN.md),
[full result](../dev/dosage_response_null_20261111/summary_23225812/results/REPORT.json),
[strata and permutation-space record](../dev/dosage_response_null_20261111/run_23225809/results/STRATA.json),
[all commands/dependencies](../dev/dosage_response_null_20261111/submission.json).

## Contract frozen before reserved values

The [reserved-test protocol](../dev/dosage_response_reserved_20261111/DESIGN.md)
and its [14 transitive source/input identities](../dev/dosage_response_reserved_20261111/SNAPSHOT.json)
are frozen before opening any of the 50 reserved RNA or DNA payloads in this
workflow. The primary candidate is gene-only. The fixed four-comparator family
is reference-only, guarded naive ratio, guarded global slope and the blend.
All use the same genes, reference hierarchy, training-only support rule and
evaluation operator.

This is a same-release reserved test within development, not a new clinical
cohort or certification that the files were never read in any historical
analysis. Known former benchmark models and related PatientIDs remain excluded.
After this test the 50 units are consumed; do not retune on them and describe a
later analysis as independently unopened.

- **Estimand:** donor-weighted 1-Mb MSE for measured relative DNA averaged over
  the fixed gene loci, and paired mean gains over all four comparators. Arm
  error and continuously truth-contrast-weighted MSE are secondary.
- **Observation model:** the frozen linear gene-response predictor, with a
  reference-only action outside the observed training RNA range. Counts are
  RNA-SeQC read abundance; no raw-count confidence model is claimed.
- **Identifying assumptions:** accurate metadata correspondence and transport
  of associations from the 119 training units. Working uncertainty assumes
  independent PatientID labels or dependence captured by lineage clustering,
  conditional on the fixed training fit. Neither is certified by a p-value.
- **Competing explanation:** lineage, culture, co-occurring alterations or
  measurement processing can predict relative DNA without a causal dosage
  response. The DNA reference itself is uncertain.
- **Distinguishing evidence:** frozen prediction on independent corresponding
  DNA, appropriate biological replication and assay resolution. Cell membership
  and spatial assignments still need truth at those respective resolutions.

Every unit remains in the primary assessment. Missing DNA remains missing and
all methods share the same measured support. Empty callable output or a failed
method fails the assessment. No unit or method may disappear to improve the gate.

The original positive-mean-gain requirement remains. This test additionally
requires positive gains against all four comparators and positive lower bounds
under both paired-donor and lineage-cluster bootstrap working models. There
are 20,000 draws per model, one fixed seed, and eight Bonferroni-protected
one-sided basic bounds. These are approximate working bands, not validated
biological coverage. Their effect, resampling unit and conditioning set are
explicit. An exact prediction tie is a failed strict superiority comparison,
not a numerical error or permission to drop that comparator.

Outside training RNA support, returning the DNA-reference action does not
identify neutral DNA. Support fractions and truth-contrast-weighted errors
report the risk of discarding real signal without introducing a convenient
CNA threshold or claiming an assay detection limit.

## Actual execution and exposure chain

Both jobs request account/QOS nslab, CPU partition, one CPU, 4 GB and one hour,
with no array or new dependency. Prediction job **23226323** completes 0:0 in
seven seconds. It fits the unchanged estimator on all 119 development units,
reads only the 50 selected RNA rows, and receives query RNA metadata only.
The raw DNA matrix is not opened in that stage.

All 11 original/support fixtures pass. Another **32** lexical estimable genes
pass independent full-training SciPy least-squares checks. Parameter serialization
and reload preserve the complete fitted digest and predictions exactly. The
zero-exceedance confidence bound is independently evaluated using 80- and
120-digit Decimal calculations before reserved RNA decoding; both converge to
the same binary64 result. The conservative rounded-up reference remains below
the fixed 5% conditional tail threshold.

All gene predictions, fitted parameters, support masks, selected IDs/order,
DNA-blind regional exports and numerical checks are hashed in the prediction
lock. Historical role labels remain in the old immutable role file. New append-
only exposure journals record when their RNA and DNA values are actually read.
This preserves history while preventing the word “unopened” from remaining an
incorrect current description after evaluation.

Dependent scoring job **23226324** is pending at this checkpoint. It must
verify the prediction lock before reading the selected DNA rows. It then checks
all regional losses independently using scalar `math.fsum`, records every paired
donor score, saves resampling indices/results and applies the frozen gate.
No superiority claim or production promotion is made before that result.

[Exact source](../dev/dosage_response_reserved_20261111/evaluate.py),
[both submitted commands](../dev/dosage_response_reserved_20261111/submission.json),
[all 11 fixtures](../dev/dosage_response_reserved_20261111/predict_23226323/NUMERICAL_CHECKS.json),
[32 full-fit references and reload check](../dev/dosage_response_reserved_20261111/predict_23226323/results/COMPUTATIONAL_CHECKS.json),
[independent confidence-bound calculation](../dev/dosage_response_reserved_20261111/predict_23226323/results/NULL_BOUND_REFERENCE.json),
[RNA exposure journal](../dev/dosage_response_reserved_20261111/predict_23226323/results/EXPOSURE.jsonl),
[prediction lock](../dev/dosage_response_reserved_20261111/predict_23226323/results/PREDICTION_LOCK.json).

## Final reserved-test result: promotion rejected

The pending statement above is a historical checkpoint. Scoring job **23226324**
subsequently completes 0:0 in 11 seconds. All 50 units remain in the result;
the primary 1-Mb assessment contains 22 lineage clusters. Every regional loss
passes the independent scalar calculation. Both RNA and DNA values for these
50 units are now consumed. They are regression-only for subsequent development.

| Frozen method | Mean 1-Mb MSE | Mean arm MSE |
|---|---:|---:|
| DNA-reference action | 0.3120896294 | 0.0551859693 |
| Guarded naive RNA ratio | 0.4968336633 | 0.1229729948 |
| Guarded global response | 0.3111966135 | 0.0546076281 |
| Gene response, primary | 0.2987860480 | 0.0520666408 |
| Blended response | 0.2987860480 | 0.0520666408 |

The fitted full-training blend weight is exactly 1.0. Its predictions are
identical to gene-only predictions, so its additional weight-fitting stage has
no demonstrated benefit. The fixed comparator remains in the family.

| Comparator | Mean paired 1-Mb gain | Donor working lower bound | Lineage-cluster working lower bound |
|---|---:|---:|---:|
| Reference | 0.0133035814 | -0.0082239254 | -0.0119752456 |
| Naive | 0.1980476153 | 0.0517211254 | 0.0298082430 |
| Global | 0.0124105655 | -0.0087935936 | -0.0124278801 |
| Blend | 0 | 0 | 0 |

Gene response beats reference on 33 units, with five exact ties. Nevertheless,
the **four-comparator/eight-band gate fails**. Reference and global bounds cross
zero independently of the blend tie. Retain the measured point improvements,
reject superiority and production promotion, and do not tune on this consumed
panel. The resampling bands are working conditional approximations; passing
their numerical reconstruction does not establish biological coverage.

[Authoritative result](../dev/dosage_response_reserved_20261111/score_23226324/results/REPORT.json),
[all paired losses](../dev/dosage_response_reserved_20261111/score_23226324/results/PAIRED_LOSSES.tsv),
[DNA exposure journal](../dev/dosage_response_reserved_20261111/score_23226324/results/EXPOSURE.jsonl).

## Secondary outcomes and independent reconstruction

The predeclared contrast-weighted endpoint weights regions continuously by
the absolute measured DNA departure from reference. It is a different estimand
from ordinary MSE and supplies neither an event threshold nor a calibrated
detection limit. At 1 Mb its mean weighted MSE is **52.0406581** for gene response,
**53.2129931** for reference, **53.1942230** for global and **49.5888550** for naive.
Thus the naive ratio is better on this secondary endpoint despite its worse
ordinary MSE. At arm resolution the corresponding values are 0.1461119,
0.1533226, 0.1524011 and 0.2322185. Preserve this resolution-dependent result.

Mean RNA support across fixed genes is **62.4158%**. A mean **39.8296%** of the
measured DNA contrast mass lies outside that RNA support. The reference action
prevents unsupported extrapolation but can discard real genomic signal. Its
fallback is not a claim of neutral DNA. This test does not identify how much
of that signal is biologically predictable from RNA.

Independent reader **23226762**, nslab CPU / 1 CPU / 2 GB / 30 min, completes
0:0 in five scheduler seconds with process peak RSS **109,888 KiB**. It
reconstructs all **160,000** saved comparison/resampling draws with scalar
`math.fsum`, reproduces all eight lower bounds and verifies exact gene/blend
identity. Maximum discrepancies are at most 2.23e-16, below bounds derived
from operation count and data scale. No refit, changed seed or new criterion
is involved. Prediction/scoring process peak RSS is 291,332 / 175,152 KiB.

The immutable summary's `primary_gate_unchanged: false` field is misnamed: its
value stores the failed gate outcome. The accompanying field note records the
canonical interpretation, **gate passed: false; criteria changed: false**.
Neither the original primary result nor the executed criteria was rewritten.

[Independent summary](../dev/dosage_response_reserved_20261111/summary_23226762/results/REPORT.json),
[summary protocol](../dev/dosage_response_reserved_20261111/SUMMARY_DESIGN.md),
[field clarification](../dev/dosage_response_reserved_20261111/SUMMARY_FIELD_NOTE.json),
[new-data feasibility and remaining measurement gaps](NEXTGEN_FEASIBILITY_2026-11-11.md).

The earlier randomization p=0.01 repeats the fixed fitting workflow, including
its internal selection, but not the development history that chose that
workflow and support correction. It is a conditional diagnostic, **not a
confirmatory p-value adjusted for method development**. Its exact Monte Carlo
tail bound concerns the frozen randomization statistic, not biological clone
false splits. The reserved-test failure is not overturned by that diagnostic.

## Current completion boundary

| Batch | Completed | Running / pending | Execution errors | Dependency status |
|---|---|---|---|---|
| Conditional-null writers + reader | 4/4 ████████ | 0 / 0 | 0 | All successful |
| Frozen RNA prediction + DNA score | 2/2 ████████ | 0 / 0 | 0 | Prediction lock verified before DNA |
| Metadata feasibility attempts | 2/2 ████████ | 0 / 0 | 0 | Schema limitation retained and corrected separately |
| Independent secondary reader | 1/1 ████████ | 0 / 0 | 0 | All saved input identities verified |

Scheduler writer peak RSS is at most 517,044 KiB; process RSS is used where
scheduler sampling reports zero. All nine jobs use account/QOS nslab. There
are no remaining HaploDose jobs in this batch. Unrelated innovation and work
jobs are untouched.

The static preservation check verifies **91** protected package/native/data/test
identities, **371/377/379** earlier ledger entries against their respective
checkpoints, **385** current ledger entries, all 11 source/prediction manifests,
18 Python parses, 11 shell syntax checks and 889 local report links at that
checkpoint. No public source changes were made in this follow-up; the public
suite was not rerun. Its historical 521-Python/four-Rust validation remains
separate from the new compute-node experiments.

[Actual scheduler accounting](../dev/dosage_response_reserved_20261111/accounting_checkpoint.tsv),
[static preservation record](../dev/dosage_response_reserved_20261111/preservation_checkpoint_v1.json),
[exact static checker](../dev/dosage_response_reserved_20261111/verify_preservation.py).

Public recovery remains 2/5, continuous-state false splits remain 14/200, and
selected correlated coverage remains 54/1,000. Stronger clone, uncertainty,
clinical, absolute-copy and broad-superiority claims remain unsupported. This
test cannot replace the missing biological null panel or suitable cell/spot
membership truth. All frozen results, past decisions and original criteria are
preserved. The full HaploDose goal remains active and unachieved.
