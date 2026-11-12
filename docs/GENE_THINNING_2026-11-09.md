# Gene-native thinning: dependence counterexample and bounded recovery

The exact distribution check passes; the proposed sufficient recovery correction
fails. Production still passes2/5 original cases. A research variant using known
simulator gene dispersion passes3/5, retaining the two previous passes and
recovering20261108;20261109/110 still fail. No runtime/default change is promoted.

## Observation law and independent evidence

The estimand here is dependence between proposal and confirmation count vectors
under repeated sampling from fixed independent gene NB laws, not a biological
false-split rate. Means and sizes are fixed across repetitions. Independence
does not mean independence conditional on the observed total, nor independence
after fitting nuisance or selecting profiles from both halves.

For aggregate mean M, variance V, beta-binomial size R and share e,
Cov(P,Q)=e(1-e)[R(V-M)-M²]/(R+1). Matching R=M²/(V-M) only makes this zero.
An independent Fraction calculation for two geometric genes with means1,4 and
sizes1,1 gives M=5,V=22,R=25/17. The low-count joint table is
[[1/10,13/200],[13/200,43/1120]], determinant−27/70000. Thus the halves are
dependent despite exactly zero covariance. Gene-native NB thinning is independent
under this fixed source model. Shared expression state remains an alternative
source of dependence even with the correct conditional split law.

Direct SciPy PMF convolution and beta-binomial probabilities, checked against
exact low-count coefficients and independent64/128-point Gauss-Jacobi quadrature,
give joint-versus-product total variation0.0167754568451 for that example.
At refined cutoff512 the tail/numerical allowance gives
[0.0167754530761,0.0167754606141]. Gene-native TV is approximately1.1e−16,
upper allowance3.77e−9. An equal-probability NB-closure control has numerical TV0,
upper allowance2.44e−10. These are distribution distances with numerical/tail
budgets, not confidence intervals or specimen error rates; budgets are not formal
interval arithmetic. The joint is never normalized after truncation. Tail bounds,
grid doubling and precision-derived budgets are saved in the gate report.

Eighteen locked contexts use all three generating populations, median true
exposure, first two mapped chr12 bins, and each of the three original failures.
Low-count joint odds ratios range0.99695748–0.99986671. Absolute probabilities
are retained; these bounded contexts do not estimate global TV or scientific error.

## Intervention and unchanged recovery criteria

The bounded research hook splits the actual retained gene matrix using simulator
size10(phi0.1), original model seed/share and global-row Philox identities, then
aggregates. The gene veto receives exactly the same gene test half. Original full
counts, search, nuisance fitting, confirmation, merge, stopping and assertions
stay fixed. Unsupported allele/spatial/null configurations are refused.

| Seed | Control pass | Gene pass | Normal false calls /90, control→gene | S1 tumour groups, control→gene | Assignment accuracy, control→gene |
|---|---|---|---|---|---|
|20261031|yes|yes|5→4|2→2|0.844444→0.836032|
|20261107|yes|yes|6→6|2→2|0.764471→0.738095|
|20261108|no|yes|10→8|1→2|0.596040→0.790514|
|20261109|no|no|2→2|1→1|0.602000→0.601202|
|20261110|no|no|8→8|1→1|0.599602→0.599602|

All S2/S3 one-group checks remain satisfied. Full six-criterion values, including
tumour detection and chr12 specificity, remain in each arm's panel/summary.json.
The two earlier passing cases lose assignment accuracy. The changed split law
also changes deterministic draw realizations; one extra pass is not a causal
attribution to independence. True dispersion is unavailable to the fitted method.
Reject this as a sufficient correction and as a deployable default. Selecting
another random split to pass the fixture would not be a justified repair.

## Audit, execution and disposition

Task: `dev/gene_thinning_20261109`; locked `DESIGN.md` predates execution.
Gate23115636, control23116789, gene23116790, audit23116791 and stronger-provenance
audit23117904 all completed0:0 on CPU innovation/nslab. Gate/audits requested
1CPU/8GB/1h; arms4CPU/8GB/2h. No arrays. Accounting and actual runtime/process
memory are in accounting_final.tsv and analysis_JOB/resources.txt. Gene arm peak
process RSS1609592KiB; elapsed1:20.43. Scheduler elapsed1:23. Actual inference
backend is CPU; execution records bind installed wheel/native/source identities.

The independent audit recounts60 original criterion values, including scientific
failures, and exactly reproduces all30 control values. Across15 fitted models it
replays every saved gene draw and checks gene/bin/veto conservation, global-row,
subset, order and chunk identities. RNG replay verifies identity, not independence.
The final audit_v2 explicitly hashes every consumed split array, fitted scorer
input and original truth table; audit_v1 is retained. The law is established by
the separate distribution calculation, not RNG identity checks.

All76 protected package/native/test/README/build files match the initial dirty
baseline. Research scripts are gene_thinning_law.py, gene_thinning_hook.py,
gene_thinning_recovery.py and gene_thinning_audit.py under benchmarks/hardening.
The previous isolated wheel validation23111358 passed395 Python and4 Rust tests,
with1 GPU-only skip; the full suite was not rerun for this research-only change.
Frozen inputs, prior outputs, criteria and historical decisions are unchanged.

Computational status: the bounded law, replay and accounting checks pass.
Identification: only fixed-known-noise conditional independence is supported;
estimated dispersion, full-count nuisance leakage, adaptive selection and
transcriptional alternatives remain unresolved. Biological validation: n=0.
No new full-workflow null/power, coverage or biological validation was run here.
The5% whole-specimen upper-bound criterion and simultaneous detection-power
requirement remain unchanged and unachieved. Further compute cannot supply
missing independent biological null specimens or suitable clone truth.
