# Rank bulk loss candidates by fitted loss-state mass

## Result

Job 23185761 completed in 21 seconds on CPU node007. The locked RNA
projection reproduces the old posterior-mean windows exactly. All 12,480
new score entries agree with independent rational median calculations.
The separate DNA stage verifies the lock, reproduces all 123 prior methods
exactly, and independently verifies the new calibration threshold and
every patient recall and neutral false-positive count. The process peak
RSS was 460,228 KiB. Slurm reported 0 KiB for the batch step.

| Fixed method | Profile Pearson, 23 patients | Loss recall, 13 patients | Neutral loss FPR, 23 patients |
|---|---:|---:|---:|
| Loss-state mass rank with native profile | 0.720291 | 0.470604 | 0.032013 |
| Native posterior mean | 0.720291 | 0.500183 | 0.057043 |
| Single-reference boundary average | 0.764817 | 0.602747 | 0.043737 |
| Segmented pseudobulk | 0.730261 | 0.399679 | 0.041244 |

Paired recall increases by 0.070925 against segmented pseudobulk and
decreases by 0.029579 and 0.132143 against the native posterior and boundary
average. Paired neutral FPR changes are −0.009231, −0.025030, and −0.011724,
respectively. Paired precision changes are +0.045057 among 11 calling
patients against pseudobulk, +0.062514 among 12 against the native posterior
and −0.014824 among 13 against the boundary average. These precision
populations differ because the methods make different numbers of calls.

The prespecified joint gate fails. The result establishes a descriptive
recall and false-positive tradeoff at each method's frozen calibration
rule. A ranking comparison at a common achieved false-positive rate
would be a different estimand. Stop this candidate's uncertainty expansion
and retain the recorded profile and call outputs. Biological replication
added is zero. Production code, frozen artifacts and scientific criteria
were preserved.

[Complete report](../dev/bulk_loss_state_rank_20261110/study_v2/score/report.json),
[all 124 patient-score methods](../dev/bulk_loss_state_rank_20261110/study_v2/score/patient_scores.tsv),
[RNA lock](../dev/bulk_loss_state_rank_20261110/study_v2/rna/RNA_LOCKED.json),
[scoring manifest](../dev/bulk_loss_state_rank_20261110/study_v2/score/run_manifest.json),
[job accounting](../dev/bulk_loss_state_rank_20261110/accounting_complete.tsv).

## Original submission and retry

The initial job 23185734 stopped in 36 seconds before computing RNA windows.
The scorer's symbol helper returns a NumPy array. The new script called a
pandas index method on that array. Job 23185761 retries after an explicit
index conversion, with a new immutable source and output directory and the
same nslab CPU allocation. The failed job, traceback and source are retained.

[Failed traceback](../dev/bulk_loss_state_rank_20261110/analysis_23185734/output.txt),
[one-file snapshot correction](../dev/bulk_loss_state_rank_20261110/snapshot_v2.json).

Job 23185734 was submitted on the CPU partition with nslab, 1 CPU, 8 GB,
and a 30-minute limit. The immutable snapshot preserves all 907 parent files.
Driver and launcher syntax checks passed before execution.

Two previously traced COAD/READ loss windows have fitted loss-state
amplitudes too shallow to reach the global cutoff on posterior mean RNA
effect. This experiment keeps the existing relative RNA profile and tests
the negative window median of the existing per-bin loss-state probability
as a separate loss-ranking score. The fitted probability describes the
working RNA model. The DNA loss labels come from the frozen DNA scorer.

Retain all 48 patients, the original 24 calibration and 24 evaluation roles,
common window support, and the patient-equal 5% neutral calibration rule.
Lock all RNA scores before opening DNA outcomes. Independently reproduce
the window medians and calibration quantile. Require all 123 old method
scores to reproduce exactly. Advancement requires improved mean paired
recall and no higher mean neutral false-positive rate against the three
fixed controls in the design. All comparator and precision outcomes remain
in the report. This is an exposed regression, with no change to the
whole-specimen clone error criterion.

[Frozen design](../dev/bulk_loss_state_rank_20261110/DESIGN.md),
[job command](../dev/bulk_loss_state_rank_20261110/jobs.json),
[snapshot](../dev/bulk_loss_state_rank_20261110/snapshot.json),
[driver](../benchmarks/hardening/bulk_loss_state_rank.py).
