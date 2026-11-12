# DNA-response development after the scientific review

The scientific goal is to improve regional prediction of independently measured
relative DNA from RNA while retaining the distinction between a predictive
association, genomic identification and validated clone membership. This is one
fixed development experiment motivated by findings 1, 3, 5, 12, 13, 14 and 17 of
the [scientific review](scientific_code_review_2026-11-10.md). Existing recovery,
false-split, power and uncertainty requirements are unchanged.

## What is implemented and what remains a scientific limitation

The [review response](SCIENTIFIC_REVIEW_RESPONSE_2026-11-10.md) records the direct
package corrections: posterior-region accumulation and complexity, uninformative
gene-veto status, residual depth balance, section-specific spatial graphs,
similarity-tree missingness, and output interpretation. Those corrections are
associated with the existing 521-Python/four-Rust compute validation and exact
one-CPU recovery replay. The public suite is not newly executed here. Protected
source identities are checked again below; successful historical execution does
not establish biological validity.

The current experiment adds no production dependency or default. It predicts a
processed relative-DNA measurement instead of treating better RNA fit as evidence
of genomic improvement. It does not provide the missing biological single-CNA
controls, same-cell membership truth, spatial registration or calibrated selected
intervals. The original production recovery remains 2/5, continuous-state false
splits remain 14/200, and correlated selected-interval coverage remains 54/1,000.

## Measurement chain and roles fixed before values

Condition/coordinate job **23225093** completes 0:0 on nslab CPU, using one CPU,
4 GB and a one-hour limit. Process elapsed time is 4.79 seconds and peak RSS is
160,960 KiB. The metadata-only gate checks model ownership, matching RNA/DNA
ModelCondition IDs and unique condition-level default profiles. It admits 171
models and excludes all 11 condition mismatches. These checks establish metadata
correspondence, not same-aliquot or same-cell correspondence.

Lexical first ModelID within each PatientID gives 169 primary patient-labelled
units; two related models remain outside the primary study. A fixed hash order
within each lineage assigns cyclic folds 0–3. Fold 0 has 50 reserved units;
folds 1–3 have 119 development units. Singleton lineages are reserved automatically,
not chosen after outcomes. PatientID labels are not certified independent donors.
This is a same-source reserved panel, not an independent clinical cohort.

Possible aliases and related individuals are excluded conservatively. In
particular, the Hs 936.T(C1) parent, LS174T parent and SCC90/same-person records
already fall in the previous exposure exclusions. Possible SCC47 and 93VU
correspondences are exclusions only, not asserted resolved identities.

Exact Ensembl IDs map to unique autosomal protein-coding GRCh38 loci in Ensembl
112. Assembly, chromosome and positive coordinates are checked. Fifteen genes
with changed biotype are excluded, leaving **17,751** fixed loci with zero-based,
half-open coordinates. This does not establish processed DNA sensitivity or
base-pair callability.

Evidence: [locked design](../dev/dosage_response_admission_20261111/DESIGN.md),
[source identities](../dev/dosage_response_admission_20261111/SNAPSHOT.json),
[complete result](../dev/dosage_response_admission_20261111/run_23225093/results/REPORT.json),
[primary roles](../dev/dosage_response_admission_20261111/run_23225093/results/primary_roles.tsv),
[gene loci](../dev/dosage_response_admission_20261111/run_23225093/results/gene_coordinates.tsv),
[alias evidence](../dev/dosage_response_admission_20261111/aliases.json),
[submission](../dev/dosage_response_admission_20261111/submission.json).

## Scientific contract for the development predictor

- **Estimand:** donor-weighted mean squared prediction error for measured linear
  relative DNA, averaged at the fixed gene loci within 1-Mb windows. Arm means
  are secondary. This is not complete physical coverage of a window.
- **Observation model:** linear prediction from per-sample RNA frequencies,
  using training lineage/library references and within-stratum response slopes.
  Counts are untransformed RNA-SeQC read abundance, not assumed UMI observations.
  No count-model uncertainty is attached to this regression.
- **Identifying assumptions:** held-out patient-labelled units do not leak into
  reference construction, coefficient fitting or weight selection; RNA/DNA
  processing and associations transport across the development splits. These
  assumptions support prediction of the processed reference only.
- **Strongest biological alternative:** regional transcription, culture effects
  and co-occurring alterations can predict DNA without being a causal dosage
  response. The processed DNA reference may itself contain systematic error.
- **Distinguishing evidence:** a frozen predictor evaluated against independently
  processed corresponding DNA, with biological replication and assay resolution
  appropriate to the claimed output. Same-cell truth remains necessary for
  stronger membership claims.

The smallest decisive analysis compares three fixed simple baselines with one
candidate on identical genes, references, support and averaging: training DNA
reference mean; naive RNA/reference ratio times that DNA mean; one pooled response
slope; and an inner-cross-validated convex blend of pooled and gene-specific
slopes. No new segmentation, smoothing, rank, start or penalty sweep is added.

The outer folds are the already locked development folds 1–3. Inner donor hash
folds select one analytic blend weight. Every reference and coefficient is fit
inside each training split. Zero predictor variance contributes exactly zero;
nonnegative projection cannot increase squared error for the nonnegative target.
Missing DNA stays missing. Infinities are invalid inputs. Unsupported training
genes and stratum/lineage/global fallbacks are recorded.

The predictor receives query RNA and RNA metadata only. It cannot use query DNA
values, callability or DNA assay metadata. Gene predictions are locked before
truth support is applied. Every method is scored with the same finite-DNA gene
operator; separate DNA-blind exports use every fixed locus. All 119 development
units remain in the primary analysis. Empty truth or predictions fail the study.

**Prespecified falsification:** reject the candidate if its primary mean error
does not beat all three fixed baselines. Differences overlapping numerical error
remain unresolved. Passing this development continuation gate would justify only
a separately frozen full-workflow null experiment, not reserved evaluation,
production promotion, genomic identification or biological error control.

## Execution, independent checks and all attempts

[Frozen protocol](../dev/dosage_response_predictive_20261111/DESIGN.md),
[eight-file source lock](../dev/dosage_response_predictive_20261111/SNAPSHOT.json),
[model](../dev/dosage_response_predictive_20261111/model.py),
[development-only loader](../dev/dosage_response_predictive_20261111/prepare.py),
[nested study](../dev/dosage_response_predictive_20261111/study.py),
[independent reader](../dev/dosage_response_predictive_20261111/audit.py).

Before freezing, static inspection corrected an exact blend-fixture expectation
from 1/2 to 1/4, rejected missing DNA infinities explicitly, and replaced a Python
3.11-only starred subscript with equivalent portable tuple concatenation. The
first login static command found no `python` executable; the `python3` parser
then exposed that syntax incompatibility. No analysis or test ran on the login
node. All five final Python sources parse and the shell launcher passes `bash -n`.

New job **23225505** uses account/QOS nslab, CPU partition, one CPU, 8 GB and a
two-hour limit, with no array. Its stages stop on failure: independent fixtures,
development measurement qualification, nested study, independent reader. It uses
the existing isolated Conda environment and native runtime from a neutral working
directory, with one numerical thread and no editable package imports. The loader
hashes original bytes but decodes numerical payloads only for the selected
development rows. Reserved and excluded payloads are not decoded.

Independent checks include SciPy design-matrix least squares, exact rational
blend optimization, feature permutation identity, constant predictors, unused
DNA metadata perturbation, a reserved-row payload that cannot parse numerically,
and outer-fold held-out DNA/callability replacement. The reader independently
checks 32 lexical estimable genes per fold with SciPy and reconstructs all
common-support regional losses with scalar `math.fsum`.

Numerical allowances derive from floating-point operation counts, design
conditioning and nonnegative accumulation bounds, with discrepancies retained.
They are computational envelopes, not biological effect thresholds. The primary
continuation differences must be separated from score-roundoff bounds; otherwise
the scientific decision stops for higher-precision investigation.

**Status at initial documentation:** job submitted; results and independent
checks are pending. No performance improvement or successful validation is
claimed at this checkpoint. [Exact submission](../dev/dosage_response_predictive_20261111/submission.json).

## Completed first experiment: rejection and traced cause

The qualification stage passes for all 119 development units. All 2,112,369 RNA
entries are finite, nonnegative integers; this does not establish raw provenance
by itself. All selected DNA entries are finite and positive, ranging from
5.847037924724008e-10 to 132.3426258799078. The reference remains processed
relative DNA, not a certified absolute copy measurement.

The writer completes all nine fixtures, qualification and the nested study.
Its original independent reader has a direct resource defect: `NpzFile` loads
and decompresses an entire array on each scalar access. Job 23225505 is stopped
after 5m26s, with its completed predictions, losses and logs preserved. The
correction loads the five arrays once, approximately 84.5 MB total, and retains
every independent formula, tolerance and decision. No completed analysis is
rerun or altered. Reader job **23225596** completes 0:0 in 18 seconds, with
245,880 KiB process peak RSS. Its new result is linked from the previously absent
expected reader path; existing outputs are not overwritten.

All **96** prespecified independent SciPy coefficient checks pass; maximum
absolute discrepancy is 1.0894063429134349e-15. Scalar reconstruction reproduces
the regional losses within their operation-derived roundoff envelopes. The
candidate's inferiority is decisively outside those numerical bounds.

| Fixed method | Mean donor 1-Mb MSE | Mean donor arm MSE |
|---|---:|---:|
| Training-stratum DNA reference | 0.295591624 | 0.0699805584 |
| Naive RNA/reference ratio | 367,077.076 | 3,033.18146 |
| One global response slope | 1,246.39900 | 9.64319633 |
| Gene-specific slopes | 9,372.40366 | 33.1732150 |
| Nested blended response candidate | 2,835.55362 | 11.2367471 |

**Decision: reject.** No reserved evaluation, production change or performance
claim. These errors measure prediction of the processed DNA reference, not
genomic clone accuracy. The no-query-RNA reference is mandatory and wins.

Cached diagnostic **23225544** then completes 0:0 in four seconds, with
313,476 KiB process peak RSS. The pending dependency was changed from the stopped
reader to its successful replacement; unrelated jobs are untouched. It
reconstructs all original predictions exactly and identifies **1,041,865 of
2,112,369** held-out RNA gene entries outside the selected training reference's
observed range; 460,650 entries have constant training RNA.

The largest donor/window error is ACH-000143, chr12:85, at five fixed gene loci.
For NTS, the query frequency is 0.00864833014 versus a constant selected training
frequency of 9.54010684e-9. The RNA-response action predicts relative DNA
109,336.6768; the processed DNA measurement is 1.01222148. Exact source columns,
counts, genomic coordinates and fitted values are retained in the diagnostic.
This traces extrapolation failure; it does not identify the biological reason
for the RNA/DNA discrepancy or certify the DNA measurement.

Replacing only out-of-range RNA actions with the existing reference-only action,
while keeping all original fitted weights, reduces the cached mean MSE to
**0.292492174**. This counterfactual supports the extrapolation explanation but
is not a new candidate result. A proper correction must refit the support rule
and weight inside every nested fold and guard the RNA baselines identically.

Evidence: [qualification](../dev/dosage_response_predictive_20261111/run_23225505/prepared/QUALIFICATION.json),
[all nine fixtures](../dev/dosage_response_predictive_20261111/run_23225505/NUMERICAL_CHECKS.json),
[original unchanged study](../dev/dosage_response_predictive_20261111/run_23225505/study/REPORT.json),
[every donor loss](../dev/dosage_response_predictive_20261111/run_23225505/study/DONOR_LOSSES.tsv),
[independent reader](../dev/dosage_response_predictive_20261111/reader_23225596/INDEPENDENT_READER.json),
[resource correction](../dev/dosage_response_predictive_20261111/READER_CORRECTION.md),
[corrected reader identities](../dev/dosage_response_predictive_20261111/READER_SNAPSHOT_v2.json),
[unchanged diagnostic protocol](../dev/dosage_response_predictive_20261111/DIAGNOSIS_DESIGN.md),
[all cause traces and counterfactual](../dev/dosage_response_predictive_20261111/diagnosis_23225544/REPORT.json),
[dependency update](../dev/dosage_response_predictive_20261111/diagnosis_dependency_update.json).

## Justified support correction, with the original criteria

The [new fixed protocol](../dev/dosage_response_support_20261111/DESIGN.md) stores
gene-level RNA min/max inside each training reference hierarchy. Outside that
observed range, each RNA-effect method uses the existing reference-only action.
It introduces no pseudocount, arbitrary cap or outcome-selected threshold. This
is abstention from an unsupported RNA effect, not DNA callability or confidence.
References, ranges, slopes and the same single convex weight are refit inside
every nested fold. All three RNA-effect baselines receive the same guard.

Job **23225716** requests one nslab CPU, 4 GB and 30 minutes. It reuses only the
qualified cached development values. All original checks and criteria remain;
two exact support fixtures are added and the independent reader is byte-identical
to the corrected first-study reader. No reserved value is decoded. Its outcome
is pending at this checkpoint.

[Frozen source](../dev/dosage_response_support_20261111/SNAPSHOT.json),
[model correction](../dev/dosage_response_support_20261111/model.py),
[nested study](../dev/dosage_response_support_20261111/study.py),
[exact submission](../dev/dosage_response_support_20261111/submission.json).

## Completed support correction and conditional validation in progress

Job **23225716** completes 0:0 in 16 seconds. All **11** fixtures pass, including
every original fixture and the two support checks. The independent reader,
unchanged from the qualified corrected reader, passes its 96 SciPy checks and
reconstructs every score. Differences exceed numerical uncertainty. Study peak
process RSS is 452,348 KiB; reader peak RSS is 253,808 KiB.

| Identically guarded method | Mean donor 1-Mb MSE | Mean donor arm MSE |
|---|---:|---:|
| Training-stratum DNA reference | 0.295591624 | 0.0699805584 |
| Naive RNA/reference ratio | 0.355782184 | 0.0925323192 |
| One global response slope | 0.294620317 | 0.0694256833 |
| Gene-specific slopes | **0.282871134** | **0.0675903490** |
| Nested blended response | 0.283507505 | 0.0676668899 |

The blended model passes its original three-baseline continuation criterion.
These are development point errors against processed DNA; no inferential
performance or biological guarantee follows. **The simpler gene-specific model
is better than blending.** Keep that result visible and do not claim that the
extra weight-fitting step adds value. This experiment is not a comparison
against the old package, cell-level pseudobulk, or the entire competitor family.
Public recovery, error/power and uncertainty outcomes remain unchanged.

[All support fixtures](../dev/dosage_response_support_20261111/run_23225716/NUMERICAL_CHECKS.json),
[every method, fold, coverage and fallback](../dev/dosage_response_support_20261111/run_23225716/study/REPORT.json),
[all donor losses](../dev/dosage_response_support_20261111/run_23225716/study/DONOR_LOSSES.tsv),
[independent numerical evidence](../dev/dosage_response_support_20261111/run_23225716/INDEPENDENT_READER.json).

The [separately frozen conditional-null protocol](../dev/dosage_response_null_20261111/DESIGN.md)
now repeats the complete fixed fitting workflow under 99 RNA/DNA pair
permutations within lineage, RNA-library and DNA-assay strata. It retains whole
RNA vectors, complete DNA structure, fixed donor folds and the same baselines.
Every draw refits references, support, response slopes and inner weights, then
passes an independent numerical reader. Three non-array nslab writers,
**23225809–23225811**, each run 33 fixed draws on one CPU, 4 GB and a two-hour
limit. Dependent summary **23225812** requests one CPU, 2 GB and 30 minutes.
No seed replacement or horizon extension is permitted.

The fixed two-target family contains the original blended model and the simpler
gene-only model. Max-statistic comparison and simultaneous exact Monte Carlo
tail bounds protect this family. Roundoff-ambiguous comparisons count as null
exceedances. All 99 draws must complete. A null draw's own continuation decision
is not an endpoint; the null reader preserves unresolved ties while retaining
every coefficient and score assertion. The observed continuation criterion is
unchanged. Its numerical qualification has already passed.

This tests conditional RNA/DNA association for a **fixed, exposed development
workflow**. It does not replay the historical choice of method, establish
exchangeability, estimate a biological false-split rate, or supply independent
clinical validation. Its Monte Carlo draws are not biological specimens. The
50 reserved units remain unopened, and neither a favourable p-value nor more
compute authorises their use for parameter selection.

[Exact null source identities](../dev/dosage_response_null_20261111/SNAPSHOT.json),
[executed workflow](../dev/dosage_response_null_20261111/workflow.py),
[all four submitted commands and dependencies](../dev/dosage_response_null_20261111/submission.json).

## Resource and preservation checkpoint

| Job | Allocation | Status at checkpoint | Recorded process peak RSS |
|---|---|---|---:|
| 23225093 | 1 CPU / 4 GB / 1 h | Complete, 0:0, 5 s | 160,960 KiB |
| 23225505 | 1 CPU / 8 GB / 2 h | Cancelled at 5m26s after writer; repeated archive decompression | 375,504 KiB completed writer; 158,316 KiB scheduler reader sample |
| 23225596 | 1 CPU / 4 GB / 30 min | Corrected reader complete, 0:0, 18 s | 245,880 KiB |
| 23225544 | 1 CPU / 4 GB / 30 min | Cause trace complete, 0:0, 4 s | 313,476 KiB |
| 23225716 | 1 CPU / 4 GB / 30 min | Support correction complete, 0:0, 16 s | 452,348 KiB study / 253,808 KiB reader |
| 23225809–23225811 | Each 1 CPU / 4 GB / 2 h | Three conditional-null writers running | Pending final accounting |
| 23225812 | 1 CPU / 2 GB / 30 min | Pending successful writers | Not allocated |

All use account/QOS nslab and the CPU partition. The cancelled job is retained
as an operational failure, not presented as a successful full execution. The
completed numerical stages and their separate successful reader are explicit.
Unrelated innovation jobs and all pre-existing running jobs are untouched.

The checkpoint verifies **91** protected package/native/data/test identities,
all **371** earlier ledger entries, six source snapshots with **38** transitive
entries, **13** Python parses, **6** shell syntax checks and **845** local links
across the four updated reports. The ledger has 377 entries; every earlier value
is unchanged. The support reader is byte-identical to the qualified corrected
reader. All 177 installed Conda manifests are inventoried; this is installed
metadata, not a newly executed environment test. The public suite is not rerun
because its protected source remains exact against the historical 521-Python/
four-Rust validation. All numerical studies above execute on compute nodes.

[Full preservation record](../dev/dosage_response_support_20261111/preservation_checkpoint_v1.json),
[installed environment inventory](../dev/dosage_response_support_20261111/environment_manifest.json),
[all nine scheduler records at this checkpoint](../dev/dosage_response_support_20261111/accounting_checkpoint.tsv).

## Completion boundary

Later completed checkpoint: scoring 23226324 consumes the fixed 50 reserved DNA
rows after the RNA prediction lock. Gene-response MSE is 0.298786 versus reference
0.312090, but uncertainty against strong baselines and the exact blend tie fail
the strict superiority gate. The secondary contrast-weighted 1-Mb endpoint
favors naive ratios. All outcomes are retained, no promotion follows, and all
50 units are now regression-only. Reader 23226762 independently reproduces all
160,000 resampling comparisons and eight bounds. The earlier pending/unopened
statements below retain their checkpoint meaning.
[Final test and limitations](DNA_RNA_RESERVED_2026-11-11.md).

Continuation:all99conditional-null draws and fourjobs now complete0:0. Both
fixed targets have0/99maxTexceedances,p.01and simultaneousMCtailupper.036575745.
This permits only conditional development continuation. The simpler gene-only
model is now frozen for the50same-source reserved units. Prediction23226323
completes0:0,with11fixtures/32additionalfull-fitSciPychecks and exact parameter
reload. Those50RNAprofiles are now decoded and predictions locked;DNA remains
unread until dependent23226324runs. The historical role lock is preserved and
current exposure recorded separately. These50must not later be called unopened.
[All null results and the frozen reserved-test chain](DNA_RNA_RESERVED_2026-11-11.md).

Computational correctness, statistical identifiability and biological validity
retain separate assessments. Even a numerically successful development model
would establish only conditional prediction against the processed reference.
The 5% whole-specimen false-split target, simultaneous power assessment,
resolution-matched validation and all uncertainty criteria remain unchanged.
Former evaluation datasets stay regression-only. Frozen outputs and historical
decisions are preserved. More computation cannot supply absent biological truth.
