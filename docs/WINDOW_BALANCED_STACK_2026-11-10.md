# Match the profile-fitting loss to represented physical windows

COMPLETE and rejected: numerical23171241 and panel23171242 both exit0:0.
The scientific gate fails despite successful execution. Both old profiles and
all16previous scores reproduce EXACTLY. The weight-only change improves one
locked case and worsens the other; no full regression, bulk transfer, weight
tuning, window-size sweep or package/default change follows.

| Locked line | Current gene-Laplace | Window-balanced | SM5 pseudobulk | MDL pseudobulk |
|---|---:|---:|---:|---:|
|2313287_STOMACH|.8155510222|.8113671996|.7487632409|.7668579269|
|42MGBA_CENTRAL_NERVOUS_SYSTEM|.8372884968|.8408883875|.7760888357|.8103601100|

The declared RNA loss improves in both cases (stomach .5973046924→.5963880017;
42MGBA .5553972870→.5552405318), beyond the numerical bounds. Better agreement
with this RNA objective therefore does not guarantee better DNA correlation.
It does not establish that either weighting universally reflects biology.
The independent rational optimum check differs by1.110223e-16,within its
derived1.080025e-11budget; dense-window replication and uniform reduction pass.
These are arithmetic checks, not biological error bounds.

Measured process peaks117008KB/numerics and333356KB/panel; scheduler elapsed
18s/23s. Scheduler MaxRSS reports0 for these short jobs, so use the retained
process measurements rather than interpreting zero as actual usage. The panel
manifest is complete with input hashes unchanged; outputs lock before DNA.
[Actual panel](../dev/window_balanced_stack_20261110/panel_v1/report.json),
[numerical report](../dev/window_balanced_stack_20261110/numerics_v1/report.json),
[RNA lock](../dev/window_balanced_stack_20261110/panel_v1/RNA_LOCKED.json),
[full manifest](../dev/window_balanced_stack_20261110/panel_v1/run_manifest.json),
[numerical resources](../dev/window_balanced_stack_20261110/analysis_23171241/resources.txt),
[panel resources](../dev/window_balanced_stack_20261110/analysis_23171242/resources.txt).

Earlier preparation checkpoint follows, retained with its original context.

Implemented one research change: weight the existing stack's absolute gene
prediction loss equally over annotated 1 Mb windows, then equally over genes
within each window. The current expert predictions, RNA preprocessing, count
bins, gene folds, penalties and scales remain unchanged. This tests whether
gene density overweights some intervals relative to the reported physical-window
endpoint. The alternative is that dense regions supply better information and
the current gene-weighted fit should remain preferable.

The estimand is regional relative RNA abundance in natural-log units against
the recorded reference. This is a descriptive working loss, not a raw-count
observation model or calibrated DNA inference. Coordinated transcription and
reference mismatch remain alternative biological explanations; corresponding
independent DNA at matching resolution is the distinguishing measurement.

Before any outcome, freeze the same two exposed development cases and require
improvement over current gene-Laplace, smoothed pseudobulk5 and Gaussian-MDL
pseudobulk in BOTH. Require all 16 previous scores and both current profiles
to reproduce EXACTLY. Save RNA outputs before DNA parsing, use identical
common support and scorer, and retain all unavailable/failing results.
If the gate fails, reject expansion without tuning windows, weights or cases.
If it passes, only a separately frozen full regression is permitted.

Independent numerical checks enumerate an exact rational two-expert optimum,
verify invariance to within-window target replication and exact reduction to
the original loss when windows contain one target. Certified optimization gaps
and operation/magnitude bounds resolve arithmetic decisions, without a new
scientific tolerance. No accuracy, calibration or superiority outcome is yet
claimed. Biological n=0 new specimens; both diagnostic lines were already used.

| Job | Role | Account/QOS | CPU / RAM / time | Submission outcome |
|---|---|---|---|---|
|23171241|Independent numerical checks|innovation|1 /4GB /30min|accepted|
|23171242|Locked two-case panel|nslab|1 /8GB /2h|accepted; afterok23171241 plus actual-report gate|

CPU v6&amd, validated old private noneditable wheel read only, captured runtime,
neutral working directory and immutable source. All 915 original source files
are exact in the new 917-file snapshot; only the research script and design are
added. No package/default/dependency/native/clone changes or frozen writes.
The existing clone jobs23167849/23169730/23170760 continue unchanged.

An initial file-management command failed because `python` is not an installed
alias; it executed no Python. The identical static/copy work succeeded using
`/usr/bin/python3`. This is recorded separately from compute execution failures.

[Frozen design](../dev/window_balanced_stack_20261110/DESIGN.md),
[implementation](../benchmarks/hardening/window_balanced_stack.py),
[source preservation](../dev/window_balanced_stack_20261110/snapshot_v1.json),
[exact job commands](../dev/window_balanced_stack_20261110/jobs.json),
[preparation attempt](../dev/window_balanced_stack_20261110/PREPARATION.md).
