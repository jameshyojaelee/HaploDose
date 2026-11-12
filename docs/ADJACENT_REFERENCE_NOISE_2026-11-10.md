# Adjacent-reference noise estimation: full recovery and a fresh error test

[Observed pattern, not a guarantee] The implemented adjacent-depth reference
estimator preserves ALL5 original recoveries, including the3 production failures,
but fails fresh error/power validation:3/10 state-null false splits versus
production0/10; exact two-clone recovery9/10 versus10/10. Reject adoption and
the final100 continuation. The original control reproduces every saved criterion,
array/table and count fold exactly. All10 CPU jobs complete0:0. Production,
original tests, thresholds, scientific criteria and frozen artifacts stay.

## A concrete noise correction, with its remaining assumptions

The preceding [depth-conditioned predictive guard](DEPTH_STATE_PREDICTIVE_2026-11-10.md)
rejects all4 known false splits but still loses original109. Do not tune its
threshold. Here change ONLY the reference noise estimator: gene rates vary
with depth within the old five-quintile reference strata. Compare nearby
reference depths to reduce that component of the factorial moment.

Within each original dataset/ref_group, sort fit-reference cells by selected-
gene count total, then unique cell ID for deterministic ties. Pair adjacent
cells disjointly; record the odd final cell as unused. No bandwidth, nearest-k,
noise floor, cap, simulator parameter or recovery-selected reference is added.
For each pair i,j and gene, accumulate

    A = (Lj/Li) Yi(Yi−1) + (Li/Lj) Yj(Yj−1) −2YiYj
    B = 2YiYj
    phi_hat = sum A / sum B; size_hat = 1 / phi_hat

For independent NB counts with known exposures Li,Lj, common phi and equal
gene rates a, E[A]=phi E[B]. If rates ai,aj differ, the exact discrepancy is
(1+phi)LiLj(ai−aj)^2. This makes the residual rate-mismatch bias explicit.
Adjacent-depth pairing narrows the depth gap without proving equal cell state.
The finite-sample ratio is not unbiased. Actual library totals are estimated
from these same counts and pairing is data-dependent, so the ideal known-
exposure identity is only a working approximation in the real pipeline.

Estimand: common residual RNA count dispersion under this reference model,
then recoverable regional RNA-profile populations. It does NOT identify pure
technical noise, DNA copy number or genomic clones. Correlated reference state,
heterogeneous gene noise, mean mismatch and estimated exposures remain competing
explanations. Independent technical measurements could separate sampling noise
from state; independently corresponding DNA, same-cell for membership, and
biological replication are necessary for genomic validation. Biological n=0.

Only declared fit-reference cells enter. No query cell, holdout reference,
confirmation array, clone label or true simulator phi chooses the estimate.
Each record saves original training row IDs, cell IDs, canonical genes, pair
indices/IDs, unmatched cells, moments and scalar size. Reference-to-feature bin
identity and observed library totals are asserted before fitting.

Arithmetic uncertainty derives from actual operation counts and moment
magnitudes, propagated through the ratio and reciprocal. It is not an SE or
confidence interval. Cross-count signal and dispersion must be positive beyond
that budget; unresolved/nonpositive values explicitly reject, without clipping.

[Prespecified design](../dev/adjacent_reference_noise_20261110/DESIGN.md),
[implemented estimator and full recovery driver](../benchmarks/hardening/adjacent_reference_noise.py).

## Independent computation and exact unchanged control

A directly calculable fixture has A=76,B=144,phi=19/36; computation agrees
exactly, within its2.19523e-13 arithmetic phi budget. Independent Fraction
expectations give equal-rate A=6/5,B=12 and changing-rate A=69/20,B=18, with
exact excess33/20 from the stated rate-mismatch formula. Row/gene permutations,
canonical tied-depth pairing, explicit odd-cell exclusion, insufficient groups,
zero cross-count support, wrong totals, duplicate IDs and nonpositive
dispersion reject correctly.
[Numerical report](../dev/adjacent_reference_noise_20261110/numerics_v2/report.json).

Before the candidate runs, the original depth_gene_score workflow with the OLD
estimator is replayed in captured runtime003. All5 original assertions and
30 criterion values,10 arrays,70 tables and15 complete count folds reproduce
EXACTLY. No tolerance or test change. The candidate starts only after this
control and numerical checks succeed.
[Strict control](../dev/adjacent_reference_noise_20261110/strict_control_v2/report.json).

## All original recovery cases, including denominator changes

| Seed | Fitted phi | Fitted size | Normal FPR /90 | Tumour TPR /510 | S1 groups | Original assignment accuracy /denominator | Every original assertion |
|---|---:|---:|---:|---:|---:|---:|---|
|20261031|0.115516819|8.656748107|0.044444|0.968627|2|0.838057 /494|yes|
|20261107|0.098852686|10.116062997|0.066667|0.988235|2|0.738095 /504|yes|
|20261108|0.103873091|9.627132394|0.088889|0.992157|2|0.762846 /506|yes|
|20261109|0.101242629|9.877262270|0.011111|0.970588|2|0.763636 /495|yes|
|20261110|0.098714674|10.130206224|0.066667|0.980392|2|0.808000 /500|yes|

Every case retains one S2 and S3 group and chr12-only reported differences.
Assignment is the ORIGINAL tumour-call-conditional criterion, interpreted with
the separate exactly-two-group/chromosome requirements, not relabeled all-cell
accuracy. The new estimator uses150 reference pairs per case. The fitted phi
values are measured reference moments, never chosen to equal simulator0.1.
For original109, normal false calls are1/90 and tumour calls495/510; membership
.763636 differs from the prior candidate's.800000, so5/5 is not a universal
membership improvement.
[All30 criterion values](../dev/adjacent_reference_noise_20261110/candidate_v2/report.json).

Posterior.95, min5 bins/min30 cells, informative count10, min8 genes, sign.01,
reference-only fitting, depth-conditioned confirmation, weighted merging and
original stopping/refinement remain. The guard rejected in the preceding report
is NOT silently added to this candidate. Existing production remains2/5.

## Fresh whole-workflow error and power

Freeze candidate/transitive hashes from source_v2, unchanged in the fresh
driver snapshot source_v4. source_v1/v3 are unexecuted drafts, not hidden
outcome-selected models. [Fresh phase7 design](../dev/adjacent_reference_noise_20261110/PILOT_DESIGN.md)
locks10 independent base seeds, matched/state/two-clone conditions and complete
production/new-estimator workflows on identical raw inputs:30 datasets,60 fits.
Reference fitting, initialization, every proposal, confirmation, gene checks,
refinement, merging and stopping rerun; truth metadata is stripped before fit.

The matched null has one CNA population; the state null adds continuous focal
chr12 transcription without a second genotype. Power adds40% second population
at delta=.6log1.5 on the SAME state background. Each condition has10 independent
replicates; paired conditions/arms are not extra independent specimens. All-
cell membership denominator600 includes incorrect/unassigned calls.

Primary candidate family:2 nulls x2 looks(pilot10/final100), exact one-sided
alpha.0125 each for simultaneous95%. Target upper error bound<=.05. At least
3/10 errors rejects excess by exact p<=1/80; any1/10 makes the best possible
final100 upper>.05, so stop that certification continuation for futility. Zero
in BOTH permits separately frozen continuation, retaining pilot. Pilot alone
cannot certify5%. All power runs finish regardless of error; power intervals
are ordinary exact two-sided95%, descriptive, without a universal power cutoff.

All60 workflows complete. Matched null0/10 errors in BOTH arms; state null0/10
production versus3/10 new candidate. The candidate's exact excess-error
p=.01150355738 is below the frozen1/80=.0125 rule; simultaneous one-sided
upper=.691501805. Even with zero further errors, its final100 upper would
be.0941614633>.05. Stop the prespecified continuation, do not increase its
horizon, tune a threshold or relabel these extra RNA groups certified clones.
Zero matched-null errors have an upper.354804988 at this pilot size, not5%.

| Method | Matched errors /10 | State errors /10 | Exact two-clone recovery /10 | Mean all600-cell membership |
|---|---:|---:|---:|---:|
|Production|0|0|10|.722000|
|Adjacent-reference candidate|0|3|9|.712167|

The ordinary exact95% recovery intervals are[.691503,1] production and
[.554984,.997471] candidate, descriptive for this fixed regime. One paired
power case favors production, none candidate;9 are shared recoveries. All
power runs finish despite null failure. No workflow returns zero groups.
The missed exact-recovery case202611097009 reports3 groups with chr12-only
differences and348/600 correct memberships, so this failure is an extra
subdivision, not an intentionally suppressed true contrast.

The preceding4/10 phase6 research errors use DIFFERENT base seeds;3/10 here
cannot be called an improvement over that old candidate. This fresh comparison
contrasts production with the whole new research candidate, not a causal
reference-estimator-only error contrast. The exact original-case control above
is the matched estimator intervention. Biological n remains0.

[Complete error, power, membership and abstention](../dev/adjacent_reference_noise_20261110/pilot_summary_v3/performance.json),
[rejection and paired recovery](../dev/adjacent_reference_noise_20261110/pilot_summary_v3/report.json).

This falsifies reference-noise correction ALONE as the requested joint repair.
The separately frozen next mechanism screen composes this estimator with the
unchanged depth-conditioned predictive safeguard on ALL5 new recoveries and
ALL3 new false splits. It does not reopen the failed unguarded phase7 final100
continuation or turn selected cases into independent calibration.

## Compute, reproducibility and preservation

Initial numerical/control/candidate jobs use CPU/v6&amd, nslab1CPU/4GB/1h,
innovation4CPU/16GB/3h and nslab4CPU/32GB/3h, respectively. Fresh6 workers each
4CPU/16GB/2h alternate accounts; dependent summary nslab1CPU/4GB/1h. No arrays
or interactive QOS. Captured runtime003, latest isolated read-only wheel,
neutral working directories and source/wheel/native/input hashes are recorded.
The previously validated456 Python/4 Rust suite is not rerun for these new
research modules; their numerical and full-workflow checks above are executed.

[Every job command](../dev/adjacent_reference_noise_20261110/jobs.json),
[all immutable snapshots and evidence](../dev/adjacent_reference_noise_20261110/).
No production package/native/original test/dependency/shared install/frozen
artifact change. All old claim values and historical decisions are retained.
The whole-specimen5% biological guarantee, independent membership truth and
full cross-mode SOTA objective remain unestablished.


| Job | Role | Account/QOS | Runtime | Outcome |
|---|---|---|---|---|
|23156825|numerics|nslab|00:00:17|complete0:0|
|23156826|strict_control|innovation|00:02:53|complete0:0|
|23156827|candidate|nslab|00:02:20|complete0:0|
|23156891|pilot_matched_00|innovation|00:01:23|complete0:0|
|23156892|pilot_matched_05|nslab|00:01:45|complete0:0|
|23156893|pilot_state_00|nslab|00:02:17|complete0:0|
|23156894|pilot_state_05|innovation|00:02:14|complete0:0|
|23156895|pilot_two_clone_00|innovation|00:02:49|complete0:0|
|23156896|pilot_two_clone_05|nslab|00:02:35|complete0:0|
|23156897|pilot_summary|nslab|00:00:08|complete0:0|

Peak scheduler RSS 1448704K (23156896.batch). Per-step records and commands are in
[accounting](../dev/adjacent_reference_noise_20261110/accounting_final.tsv) and
[resource summary](../dev/adjacent_reference_noise_20261110/resource_summary_final.json).
