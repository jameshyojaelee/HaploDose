# Runnable REFCON comparison for single-cell regional profiles

The public REFCON release runs on the cluster and reproduces all four
published A375 metrics. This verifies the released comparator on its own
fixture. The matched scONE analysis uses the same rounded expression and
annotated normal reference for HaploDose and REFCON's reference-refined
option. It uses one exposed specimen and published DNA depth profiles.
The result has a narrower scope than genomic membership or independent
validation.

Later development checkpoint: the original default-primal scONE fit/scorer
23188436/23188477 completes. Whole/cell-specific mean Pearson agrees with
the dual comparison at the displayed precision; the full arrays have a
maximum 1.63×10⁻⁸ discrepancy, so they are not declared byte-equivalent.
A separate anchored RNA candidate reaches 0.573712/0.297868, compared with
pseudobulk 0.561258/−0.048242 and REFCON with reference 0.615822/0.376007.
Only 161/390 cells improve whole Pearson over pseudobulk. This remains
one exposed specimen and no default is promoted.
[All methods, matched controls and default-primal verification](SCONE_REGIONAL_SHRINK_2026-11-10.md).

## A375 release reproduction

The run pins REFCON source to `ad70105cc18fa29e0d9089e7fcd6caba2405a0fb`
and checks the downloaded Zenodo archives against repository SHA-256 and
archive API size and MD5. It verifies both checkpoint identities. RNA-only
input retains 96 cells and 14,894 genes in release order. Both checkpoints
run with forward and reverse passes. Both checkpoint results and the
ensemble preserve the input cell and gene identities. All predictions
are locked before opening the released DNA array.

| Native endpoint | Independent recomputation | Published at displayed precision | Evaluable cells |
|---|---:|---:|---:|
| Pearson | 0.469635 | 0.470 | 64 |
| Spearman | 0.479025 | 0.479 | 64 |
| Loss AUROC | 0.734419 | 0.734 | 64 |
| Gain AUROC | 0.758547 | 0.759 | 64 |

The independent scorer uses tied ranks and the rank-sum AUROC formula.
It reproduces the author scorer at four decimals and the released table
at three. The 64 cells come from one cell line. The author-provided DNA
file uses positional identity. This reproduces the release score and
leaves the underlying DNA callability and correspondence requirements
in force. No HaploDose A375 comparison is made with the release's fixture,
which provides no matched normal expression pool.

Job 23188381 completed `0:0` in 3 minutes 10 seconds on `ne1dc6-010`.
The `nslab` CPU request was 4 CPUs, 24 GB and 3 hours. REFCON inference
took 87.1 seconds across both models and 439,180 KiB peak process RSS.
The scheduler recorded 501,108 KiB peak job RSS. PyTorch 2.14.0+cpu ran
with four threads in the existing read-only Conda environment.

[REFCON source](https://github.com/ciceklab/refcon/tree/ad70105cc18fa29e0d9089e7fcd6caba2405a0fb),
[release scoring protocol](https://github.com/ciceklab/refcon/blob/ad70105cc18fa29e0d9089e7fcd6caba2405a0fb/docs/reproduce.md),
[Zenodo record](https://doi.org/10.5281/zenodo.23534511),
[A375 result](../dev/refcon_regression_20261110/release_23188381/A375_REPORT.json),
[prediction lock](../dev/refcon_regression_20261110/release_23188381/A375_LOCKED.json),
[compute environment](../dev/refcon_regression_20261110/release_23188381/environment.json),
[submission](../dev/refcon_regression_20261110/release_submission.json).

## Matched scONE regression

The comparison targets all 390 annotated tumor cells with published same-cell
processed DNA, with 450 annotated normal cells available for the two
reference-based methods. The regional operator uses fixed 5-Mb windows
split at GRCh38 arm boundaries. Every method uses the same finite support
at scoring. For the cell-specific endpoint, profiles are put on a common
mean-one scale before subtracting the other tumor cells' regional mean.
The analysis keeps an unscaled diagnostic as well. The windows define
an evaluation operator and supply no independent DNA detection limit.
The release describes rounded RNA values, which remain reconstructed input
for HaploDose's descriptive APIs. The fit and score use frozen, separate
artifacts. The previous comparator family and other historical results
are preserved.

The original primal scONE score is pending. Job 23188436 fits HaploDose on the locked RNA
input. REFCON jobs 23188473–23188476 cover the four CPU shards and the
32 repeated probe cells. The two earlier pending jobs, 23188437 and
23188438, were canceled before execution when the L40S had a scheduler
start estimate almost two hours out. The replacement scorer 23188477
requires the five successful fit jobs and exact cross-shard probe parity
before opening the scONE DNA matrix. All submitted jobs use `nslab`.

The four REFCON shards completed successfully. Separate job 23188921
reassembled all 840 cells with exact agreement for the 32 repeated cells
and preserved every source axis. The original HaploDose primal fit still
runs. The completed dual-solver replay used these locked REFCON predictions.
[Exact shard check](../dev/refcon_regression_20261110/probe_23188921/SHARD_PARITY.json),
[solver and replay record](STACK_SOLVER_SCALING_2026-11-10.md).

### Completed dual-solver comparison

The five methods were evaluated on 390 cells, 572 fixed windows and 22
autosomes. The common finite gene/cell mask contains 11,260,470 pairs.
These are within-specimen descriptive means. The biological unit is one
exposed specimen, with processed same-cell DNA relative depth as the
comparison measurement.

| Method | Mean Pearson | Mean Spearman | Mean cell-specific Pearson |
|---|---:|---:|---:|
| REFCON, reference-free | 0.432251 | 0.417711 | 0.371134 |
| REFCON, reference-refined | 0.615822 | 0.557644 | 0.376007 |
| HaploDose averaged stack | 0.419387 | 0.412188 | 0.226796 |
| Naive individual RNA ratio | 0.387827 | 0.393342 | 0.198905 |
| Leave-one-cell-out pseudobulk | 0.561258 | 0.399008 | −0.048242 |

The HaploDose stack improves both Pearson endpoints over the naive
individual ratio. It trails pseudobulk for whole-profile Pearson and
both REFCON options for whole-profile and cell-specific Pearson.
Pseudobulk carries the shared profile, so its leave-one-out residual
should be interpreted in the context of the centring operator. These
correlations establish neither event accuracy nor calibrated membership.
The cohort remains regression-only.

All four CPU fits completed `0:0` in 4m53s–4m57s, with 419,168–423,172
KiB peak scheduler RSS. The repeated first-cell HaploDose probe was exact
across all shards. The scorer completed `0:0` in 43s, with 842,064 KiB
peak scheduler RSS. Prediction hashes were verified before the scoring
stage opened the DNA data.

[Five-method report](../dev/stack_scaling_20261110/scone_score_23188976/SCONE_REPORT.json),
[per-cell scores](../dev/stack_scaling_20261110/scone_score_23188976/scored_cells.tsv),
[prediction lock](../dev/stack_scaling_20261110/scone_score_23188976/ALL_PREDICTIONS_LOCKED.json),
[HaploDose shard parity](../dev/stack_scaling_20261110/scone_score_23188976/SHARD_PARITY.json),
[Slurm accounting](../dev/stack_scaling_20261110/accounting_scored.tsv).

[Frozen design](../dev/refcon_regression_20261110/MATCHED_DESIGN.md),
[CPU and scoring amendment](../dev/refcon_regression_20261110/MATCHED_AMENDMENT.md),
[initial submission log](../dev/refcon_regression_20261110/matched_submissions.jsonl),
[replacement submission log](../dev/refcon_regression_20261110/cpu_submissions.jsonl),
[accounting checkpoint](../dev/refcon_regression_20261110/accounting_checkpoint.tsv),
[scONE measurement feasibility](SCONE_FEASIBILITY_2026-11-09.md),
[DNTR correspondence](DNTR_CORRESPONDENCE_2026-11-09.md).
