# Faster masked-gene RNA profile stacking

The same conditional stacking loss can be solved using one linear-program
constraint per expert. In the fixed scONE RNA-only probe, this reduces
full fit time from 24.20 to 2.90 seconds on one CPU. Both formulations
select the same stack. The maximum profile difference is 2.22×10⁻¹⁶.
The package adds `solver="dual"` as an opt-in mode. The primal solver
remains the default to preserve historical choices when weights are tied.

## Cause and correction

Completed full-fixture follow-up: primal 23188436 finishes in 2h37m17s,
with scorer 23188477 also complete. Comparison 23190558 finds identical
geometry, both simple-baseline matrices and 1,560 non-stack cell-score
records. Stack profiles differ by at most 1.63×10⁻⁸ and cell scores by
1.24×10⁻⁹. Preserve both outputs and the opt-in/default distinction.
[Full default comparison](../dev/scone_regional_shrink_20261110/run_23190558/REPORT.json).

Profiling the lexically first eligible scONE tumor cell gives 28.08 seconds
over 29,259 genes and 572 bins. The three HiGHS calls account for 21.69
seconds. The existing primal problem has one auxiliary variable and two
inequalities per masked gene. The dual has one box-constrained variable per
gene and one inequality per expert. The solver's inequality multipliers
recover simplex weights. The existing independent loss and feasible dual
checks still apply. An unsuccessful dual solve falls back to the primal.

The exact small two-expert counterexample passes. The real probe then
reproduces the saved primal output exactly and checks the two formulations
on the same inputs, CPU, SciPy version and thread count. All three primal
losses agree exactly at binary64 precision in that probe.

| Experts | Eligible genes | Primal solve | Dual solve | Maximum weight difference |
|---:|---:|---:|---:|---:|
| 6 | 29,257 | 6.64 s | 0.147 s | 3.89×10⁻¹⁶ |
| 7 | 29,257 | 7.83 s | 0.166 s | 5.86×10⁻¹⁶ |
| 8 | 29,257 | 7.31 s | 0.178 s | 5.55×10⁻¹⁷ |

Different weights can give the same objective when the masked expert
predictions are linearly dependent. Full profiles can then differ. The
observed numerical agreement is specific to the probe and supplies no
general profile-equivalence or biological claim. The optimization changes
neither the loss nor the independent numerical decision gate. The existing
scONE primal job continues with its already installed wheel.

[Profiling output](../dev/stack_scaling_20261110/profile_23188910/profile.txt),
[objective derivation and design](../dev/stack_scaling_20261110/DUAL_DESIGN.md),
[probe results](../dev/stack_scaling_20261110/dual_23188920/report.json),
[public solver](../python/haplodose/stacking.py).

## Compute and preservation

Profiling job 23188910 completed `0:0` in 51 seconds with 352,520 KiB peak
scheduler RSS. The same-node solver probe 23188920 completed `0:0` in
31 seconds, with 279,924 KiB peak process RSS. The separate REFCON shard
check 23188921 completed `0:0` and found exact prediction agreement for
all 32 repeated cells across the four shards. These jobs use `nslab` CPU.
No DNA was opened by these jobs.

Environment 23188947 completed in 3 minutes 47 seconds. Wheel and suite
job 23188948 completed `0:0` in 2 minutes 56 seconds. The installed wheel
passed 519 Python tests and four Rust tests, with one GPU-only skip and 15
HiGHS option-forwarding warnings. Installed replay 23188949 requires
exact old profiles and diagnostics on the fixed scONE probe and one
prespecified case from each of the four archived pooled-profile chunks.
All five exact primal checks passed. The dual matches the research probe
exactly. On the four pooled cases, selection decisions and objective
values agree, the largest profile difference is 9.99×10⁻¹⁶, and the
full fit is 2.2–2.8× faster. The installed scONE probe runs in 2.90 s
versus 24.15 s with the primal. Replay job 23188949 completed `0:0` in
1 minute 19 seconds with 274,848 KiB peak process RSS.

The separate scONE dual run completed four fit jobs, 23188972–23188975,
and its dependent scorer, 23188976. It uses the same locked RNA geometry
and processed-DNA scorer as the original regression. The shared probe
matches exactly across shards. The 390-cell score trails pseudobulk on
whole-profile Pearson and reference-refined REFCON on both Pearson
endpoints. The original primal fit still runs.
[Five-method scores](REFCON_REGRESSION_2026-11-10.md#completed-dual-solver-comparison).

[Source snapshot](../dev/stack_scaling_20261110/snapshot.json),
[installed validation submissions](../dev/stack_scaling_20261110/installed_submissions.jsonl),
[installed wheel](../dev/stack_scaling_20261110/validation_23188948/installed_artifact.json),
[Python tests](../dev/stack_scaling_20261110/validation_23188948/pytest.txt),
[Rust tests](../dev/stack_scaling_20261110/validation_23188948/cargo_test.txt),
[installed replay](../dev/stack_scaling_20261110/installed_23188949/report.json),
[scONE solver replay design](../dev/stack_scaling_20261110/FAST_SCONE_DESIGN.md),
[scONE replay submissions](../dev/stack_scaling_20261110/fast_scone_submissions.jsonl),
[REFCON shard check](../dev/refcon_regression_20261110/probe_23188921/SHARD_PARITY.json),
[Slurm checkpoint](../dev/stack_scaling_20261110/accounting_checkpoint.tsv),
[installed and follow-up job accounting](../dev/stack_scaling_20261110/accounting_installed.tsv).
