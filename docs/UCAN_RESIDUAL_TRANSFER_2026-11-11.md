# Raw-count U-CAN transfer: failed candidate, corrected arm scoring

The scientific goal remains accurate RNA-derived regional genomic profiles and
supported CNA-defined populations, with independently established error and
uncertainty limits. The unchanged residual-precision candidate passed its fixed
simulation replication but fails this locked real-count regression:16/72 primary
comparisons and0/12 state-subtraction ablations. Reject general transfer and
default promotion. A separate, necessary benchmark correction makes arm summaries
and their tied ranks independent of reduction rounding and input dtype. It does
not change the failed primary decision, frozen outputs or scientific criteria.

## Scientific contract and fixed experiment

This addresses findings5,12–14 and17 of the
[November10 review](scientific_code_review_2026-11-10.md). Its predecessor is the
[unchanged-estimator replication](WITHIN_REGION_REPLICATION_2026-11-11.md).
The [frozen design](../dev/ucan_residual_transfer_20261111/DESIGN.md) was written
before new fitting or DNA scoring. These two already exposed diagnostic tumours
are regression-only:CRC.SW.U0001.T andCRC.SW.U0008.T. Biological query n=2;
120 distinct normal donors do not turn this into122 independent tumours.

The estimand is a centered regional natural-log RNA response relative to the
specified fit normals, after subtracting state predicted from within-region gene
contrasts. Counts are raw featureCounts paired fragments, not reconstructed
TCGA expression. Region-excluded residual variance estimates working precision;
it does not establish sampling variance, confidence coverage or genomic identity.

The model assumes informative within-region state markers, sufficiently uniform
regional gene dosage response, and transport of the normal state-to-region
relationship to tumour. The strongest competing explanations are tumour-specific
transcription, cell composition, reference mismatch and gene-specific dosage
response. Independent DNA on corresponding material distinguishes profile
compatibility; it cannot identify which of those mechanisms caused a correction
failure. FACETS segment boundaries, purity and ploidy are uncertain measurements.
These comparisons establish neither absolute copy number nor cellular membership.

All original RNA-only gene masks, ordered genes,1412 bins and genomic blocks were
retained. Lexically ordered normal metadata was permuted with seed2026111107:
60 fit donors and60 selection donors, with no all-normal refit. The old RNA mask
used all120 normals; selection loss is conditional on that shared preprocessing.
Molecule halves use seed2026111108 and conserve every input count exactly.
Fit/selection roles, identifiers and order are saved in
[the normal-role table](../dev/ucan_residual_transfer_20261111/run_23206386/results/normal_roles.tsv).

The frozen method uses normal-only ridge prediction, a fixed14-point path and
region-excluded normal residual weights. No query RNA or DNA selects the path,
weights, donors or penalty. A driver-local feature cache must reproduce all28
archived profile arrays and14 parameter dictionaries exactly before fitting.
That check and inherited independent numerical fixtures pass.

Both specimens' full24-method RNA families were locked before either DNA table
or old DNA score was decoded. Six old base estimators each retain raw, width5
median and the same RNA-penalized independent Gaussian partition. All36 old
arrays replay exactly. New corrected and matched weighted-raw profiles receive
the same three operators. The full transfer gate requires strict superiority
over all six matched baselines, both primary correlations, all operators and
both specimens:72 comparisons, plus12 state-subtraction ablations. No epsilon
defines a win, and no comparator disappears. This gate remains unchanged.

## Measurement and execution chain

The source verification reconstructs each selected gene count from the original
gene-ID matrix, unique Ensembl112 annotation and both tumour-specific exon files.
Counts, samples and positions match the old canonical input exactly. Used-gene
positions are checked against one-based input and zero-based half-open output.
Every normal donor is unique. Original featureCounts/samtools runtime versions
remain missing; successful reconstruction does not fill that provenance gap.

The four linked records are:

1. [Input identities and source hashes](../dev/ucan_residual_transfer_20261111/run_23206386/results/INPUTS_LOCKED.json),
   [gene mapping](../dev/ucan_residual_transfer_20261111/run_23206386/results/input_identity_map.tsv)
   and [selected genes](../dev/ucan_residual_transfer_20261111/run_23206386/results/genes.tsv).
2. [Executed source manifest](../dev/ucan_residual_transfer_20261111/source.sha256),
   [unchanged estimator](../dev/ucan_residual_transfer_20261111/model/conditional.py),
   [driver](../dev/ucan_residual_transfer_20261111/study.py) and
   [fit environment](../dev/ucan_residual_transfer_20261111/run_23206386/results/fit_environment.json).
3. [Candidate RNA lock](../dev/ucan_residual_transfer_20261111/run_23206386/results/CANDIDATE_RNA_LOCKED.json)
   and [both complete RNA families](../dev/ucan_residual_transfer_20261111/run_23206386/results/ALL_RNA_LOCKED.json).
4. [DNA paths, coverage and every primary decision](../dev/ucan_residual_transfer_20261111/run_23206386/results/REPORT.json),
   [all144 original score rows](../dev/ucan_residual_transfer_20261111/run_23206386/results/scores.tsv),
   [all108 historical comparisons](../dev/ucan_residual_transfer_20261111/run_23206386/results/historical_scores.tsv)
   and [independent reader](../dev/ucan_residual_transfer_20261111/reader_v5.py).

DNA is reconstructed independently of RNA grouping. Whole-gene midpoint coverage
and complete-bin requirements give1410/1412 callable bins inU0001 and1409/1412
inU0008. No candidate-specific support reduction is permitted. Independent
`bisect` interval lookup and standard-library median/log calculations reproduce
every DNA profile value exactly. No inconsistent count, barcode, assembly or
coverage evidence arose that justified BAM reprocessing. Read inspection was
not run. The tumour-specific DNA material correspondence remains limited by
the existing provenance, including theU0008 TB-labelled DNA file; no new
specimen correspondence or truth certification is claimed.

## Results and mechanism

Below, every row uses the same Gaussian partition operator and callable support.
These are exposed-specimen point estimates, not patient-level superiority tests.

| Method | U0001 Pearson / Spearman | U0008 Pearson / Spearman |
|---|---:|---:|
| Naive ratios | .780057 / .766850 | .915604 / .590964 |
| Native posterior | .802236 / .795957 | .959666 / .623996 |
| Posterior-region profile | .839552 / .840990 | .919602 / .630539 |
| New residual precision, no state subtraction | .838318 / .830213 | .938366 / .644356 |
| Same candidate with state subtraction | .694174 / .634598 | .791250 / .626794 |

The full candidate wins16/72 primary comparisons and loses all12 matched
subtraction comparisons. Weighted raw profiles have promising individual point
results but do not establish superiority over the full comparator family.
Do not select a different estimator per specimen or promote the raw ablation
retrospectively. The prior successful96/96+12/12 simulation replication remains
true under its declared laws; it did not establish tumour transport.

Normal selection chooses ridge index13,relative alpha0.000244140625. Selection
loss decreases from0.01224343 without subtraction to0.00324068 at that index.
Normal RNA prediction therefore improves while tumour DNA compatibility worsens.
That is a direct counterexample to treating normal prediction improvement as
evidence of improved genomic inference.

Before partitioning, removed components have RNA-effect RMS .171607/.170564
and DNA Pearson .407448/.333657 inU0001/U0008. Raw DNA Pearson falls
.597887→.501410 and .680755→.641265. Under partitioning, the profile differences
have DNA Pearson .619901/.537279; these include the nonlinear operator response
and are not the original linear predictor itself. [Inference] This supports
loss of DNA-aligned information through subtraction. It does not establish
whether tumour state transport, gene dosage response, admixture or uncertain
DNA correspondence is responsible. Lowering the ridge penalty or using DNA to
choose normals would not resolve that identification problem.

The smallest decisive biological follow-up would measure independently matched
DNA and RNA across replicated normal/tumour material with composition metadata,
then test the frozen state predictor at the declared resolution. Without that
evidence, retain the simpler RNA path and withhold genomic interpretation of
the subtraction. More simulation draws do not supply missing transport evidence.

## Necessary correction: exact arm summaries and ranks

The first independent reader confirms all84 primary decisions using integer
cross-products, independently coded average-tie ranks and80/160-digit Decimal
correlations. All576 numerical references are retained, but secondary scores
do not all reproduce. This distinction must remain explicit.

Two mechanisms were traced:

- Pandas can retain float32 through group medians. Averaging the middle pair
  rounds to float32 before correlation, whereas exactly the same observations
  promoted to float64 retain that midpoint. The independent median-score
  discrepancy reaches2.8722794e-9. This is not explained by conditioning.
- Rounded group means can create tiny differences between mathematically
  identical arm means with different bin counts. Spearman ranks treat those
  differences as ordered observations; the largest observed score change is
  .0453170645. An arbitrary rank tolerance would redefine ties without a
  scientific justification.

The new [arm scorer](../benchmarks/hardening/arm_metrics.py) forms exact rational
means/medians of supplied binary observations and exact average-tie ranks.
Correlation displays use80-digit Decimal arithmetic. Both live
[U-CAN](../benchmarks/hardening/ucan.py) and
[posterior-transfer](../benchmarks/hardening/ucan_posterior_transfer.py) benchmark
drivers now use it for arm scores and record its transitive source hash.
Raw common-bin scoring, fitted package outputs, thresholds and frozen results
are unchanged. This correction concerns numerical summaries of the measured
values, not biological truth or the uncertainty of an arm effect.

[Ninety-six corrected secondary score rows](../dev/ucan_residual_transfer_20261111/reader_v5_23206559/results/corrected_arm_scores.tsv)
sit alongside all historical scores; nothing is overwritten. The separately
implemented exact-mean reader verifies192 mean-metric references, including
bitwise agreement of exact Spearman displays. Every primary decision remains
16/72+0/12, with maximum independent primary display discrepancy4.4408921e-16.
The final reader deliberately reports `scoring_fully_verified:false` for the
old scoring path; successful execution is not a waiver of its discrepancies.

The [eight new invariant tests](../tests/test_arm_metrics.py) pass from a frozen
snapshot on a compute node. They check independently derived35/38 tied-rank
correlation, float32 promotion of identical observations, within/between-arm
permutations, nonfinite fixed-support rejection and the live wrapper from a
neutral directory. The original tests and criteria were not edited. The prior
521-Python/four-Rust package validation is historical, not rerun here: no public
package code changed.

## All attempts, resources and preservation

All seven jobs used nslab account/QOS,CPU,one CPU. Main request:16GB/3h;
each reader/test request:8GB/1h. Every batch header was announced before
submission. No arrays, interactive QOS, GPU, shared installs or running-job
changes. [Exact commands, headers and GNU-time memory](../dev/ucan_residual_transfer_20261111/jobs_complete.json)
and [scheduler accounting](../dev/ucan_residual_transfer_20261111/accounting_complete.tsv)
retain every attempt. Short-job Slurm RSS is sometimes0; GNU-time process RSS
is used below, not interpreted as allocated memory or total aggregate usage.

| Job / stage | Progress | Outcome | Elapsed | Peak process RSS,KiB |
|---|---|---|---:|---:|
| 23206386 prepare/fit/score | 4/4 ████████ | Complete0:0; scientific gate fails | 3m51s | 403448 |
| 23206425 reader1 | 1/1 ████████ | Complete0:0; primary verified, secondary discrepancy found | 5s | 110552 |
| 23206436 reader2 | 0/1 ░░░░░░░░ | Failed1:0; unconditioned arithmetic assertion | 7s | 109036 |
| 23206446 reader3 | 0/1 ░░░░░░░░ | Failed1:0; conditioning does not explain aggregation discrepancy | 7s | 109024 |
| 23206517 reader4 | 0/1 ░░░░░░░░ | Failed1:0; retained tables, report serialization rejects NumPy bool | 7s | 109008 |
| 23206559 reader5/corrected scorer | 1/1 ████████ | Complete0:0; exact primary and corrected-mean checks | 10s | 111508 |
| 23206566 focused invariant suite | 8/8 ████████ | Complete0:0;8 tests pass | 16s | 127316 |

Running0,pending0,complete4,failed3 in this chain. Downstream fit/score stops
on executable failure. Reader failures did not overwrite the main run; all
versions, logs and the reader4 incomplete JSON remain. Reader1's broad numerical
diagnostic was not a pass criterion. Reader2 incorrectly budgeted only final
arithmetic; reader3's centering bound still failed. Reader4 then retained every
diagnostic flag instead of stopping at the first one, and reader5 repaired only
serialization, added the exact scorer and kept the old-score failure flags.
No failed assertion was relabelled a scientific pass.

One login static-check attempt used unavailable bare `python`; the isolated
Conda interpreter subsequently parsed the files successfully. No scientific
execution ran on the login node. Python/scipy/pandas versions and CPU nodes
are recorded per stage. The main fitting uses the simulation's unchanged Conda
interpreter; preparation/scoring/tests use the validated noneditable wheel
environment under its shared read lock. No environment was modified.

The [integrity record](../dev/ucan_residual_transfer_20261111/source_integrity_complete.json)
confirms all91 previously validated package/native/data/test identities and all
execution source manifests are unchanged. New benchmark/test files match their
validation snapshot. Before-edit benchmark copies and the dirty tree inventory
remain under the task directory. The claim ledger adds one entry while retaining
all280 previous values. Frozen artifacts and the original review are unchanged.

The full goal is not achieved. Original production recovery remains2/5;
continuous-state false splits14/200 fail the5% bound, and correlated selected-
interval coverage remains54/1000. This experiment supplies no new clone-recovery,
membership, spatial or uncertainty validation, no independent evaluation cohort
and no SOTA claim. The narrow implementation correction is accepted; the
scientific state-subtraction transfer is rejected.
