# Test whether regional mean normalization repairs recovery/state separation

**Complete result: 7/9, reject this correction as sufficient.** The same true
20261031 contrast remains negative and the same state control remains unresolved.
The independent reader verifies all 27 score arrays, nine decisions, 36 full-PMF
probes and all nine source/exposure reconstructions. Fifteen nslab jobs complete;
one initial exact-replay failure remains preserved. No production promotion,
scientific criterion change or frozen-result change follows. Completion details
below follow the earlier execution checkpoint, which remains historical.

The scientific goal remains independently supported genomic dosage profiles,
CNA-defined cell/spot populations and their spatial distribution, with the
original error and power criteria. This experiment changes the recovery guard's
regional observation mean. It does not change thresholds, split counts,
membership priors, biological labels or the interpretation of frozen results.

## Mechanism and scientific contract

The preceding joint-covariance candidate fails real-cell profile performance.
Earlier whole-cell confirmation and reserved allele observations also fail to
separate all state-driven distinctions. Repeating those methods or treating
reserved SNPs as DNA evidence would not supply a justified recovery correction.
The present intervention instead tests a discrepancy between the assignment
model's compositional means and the guard's unbounded exponential means.

For a baseline regional RNA share w and latent regional effect eta, replace
`E_g*exp(eta)` by

    mu_g = E_g*exp(eta)/(1-w+w*exp(eta)).

This conserves a two-block baseline total **expectation** when the outside
expectation receives the same normalizer. It bounds the expected regional
fraction. It does not magically turn a product of NB scores into a joint count
distribution conditional on its own observed total. At independently fixed
exposures it defines a gene-count NB model. The actual historical nuisance and
selection here are conditioned upon, so its result remains a conditional
mechanism comparison, not independent confirmation or calibrated specimen error.

The estimand is confirmation-RNA predictive advantage of a discrete regional
population mixture over continuous state. In both families, eta has a common
linear outside-event depth trend and shared Gaussian variance; the alternative
adds an effect separation and depth-independent prevalence. Independent gene-NB
noise, Gaussian state, transferable frozen expectations and the two-block
normalization are working assumptions. Regional transcription, gene dependence,
depth-model error and reference mismatch remain alternatives. A transcriptional
state and a CNA can still give exactly the same RNA law. No softmax resolves
that identification problem. Independent same-cell DNA and biological controls
are needed for the genomic claim; technical measurements can test noise.
Biological n=0 in this exposed nine-proposal simulation screen.

Set w=A/(A+B), where A is the frozen expected discovery event count and B is the
discovery count over genes outside every selected event bin. This is a plug-in
baseline exposure estimate, not measured absolute depth. The old binomial split,
source counts, row/gene IDs and counter domain are reproduced without any new
random allocation. All old discovery fields remain exact. Neither archived q
from a different split nor new confirmation values estimate w. All nine share
ranges are modest, approximately 0.0148–0.0306; they are recorded per case.

The normalized log response is h. Its derivatives are h'=1-q and h''=-q(1-q),
where q is the shifted regional fraction. Therefore its log-likelihood Hessian
is `h'^2*H_raw + h''*S_raw`. The variance-zero heat derivative must use this
complete Hessian and the squared transformed score. The implementation does so,
including the analytic zero-variance case. No positive variance floor is added.

The same genes in discovery and confirmation share sampling intensity. As in the
previous experiment, conditional predictive evidence is the full-count marginal
model difference minus the discovery marginal difference. The binomial allocation
factor cancels. Its gene-NB shape is not scaled by the split fraction. The same
fixed w appears in both marginals. New fits retain the two/seven starts, optimizer,
24/48/96 quadrature, q initialization and minimum component support of 30.

All five true contrasts must resolve positive with supported components, and all
four state controls must resolve negative. All nine are required. Numerical
resolution remains twice the largest consecutive-order score discrepancy plus
operation/magnitude rounding budgets; it is not a confidence interval or rigorous
global integration enclosure. An unresolved state score does not become a
successful rejection by rounding. A screen pass would only permit separately
frozen complete-workflow recovery and fresh false-split/power testing with the
unchanged 5% simultaneous specimen bound requirement. Original recovery 2/5,
continuous-state false splits 14/200 and selected coverage 54/1,000 remain open.

## Numerical prerequisite and preserved failed execution

The initial prerequisite 23207181 passes the independent normalized-probability,
five-gradient and variance-boundary checks. Rational discrepancy is 3.05e-16;
independent SciPy full-PMF/adaptive Gaussian discrepancy is 1.67e-16. Their
operation/tail budgets are approximately 1.15e-9. Finite differences, exact
two-block expected-total conservation, finite extreme effects, component/gene/
cell permutation, nested null, covariate reparameterization and the exact
zero-share exponential limit also pass.

It then fails exact archived replay before any new fit: 29/599 values in state
case index 7 differ by at most 2.88658e-15. The job is terminal FAILED1:0 after
39 seconds on ne1dc4-003, a v4 Intel node. Its archived case ran on v3 Intel
ne1dc3-022. This is an actual failed equality gate, preserved with its nine-file
source snapshot. No array or tolerance is removed.

Preparation retry 23207217 separates data creation from hardware-specific
replay. It completes in 14 seconds. Replay 23207237 runs the eight original v4
cases on ne1dc4-001; replay 23207238 runs index 7 on its original ne1dc3-022.
Both the unchanged old model and the new zero-share model must match every
archived array exactly. They do: 24+3 arrays, with no refitting. Qualification
23207239 checks the complete case/order set and succeeds before new fitting.
This supports an execution-context explanation for this particular discrepancy;
it does not establish universal cross-hardware bitwise reproducibility.

## Execution checkpoint

Nine new case jobs 23207253–23207261 and summary 23207262 are submitted after
qualification. Reader 23207272 waits for the successful complete summary. No
scientific gate is claimed from partial case execution. Every fit locks all
three orders before its confirmation counts are decoded for scoring.

All submissions use nslab account/QOS and CPU. Preparation/replays/qualification:
1 CPU/8 GB/1 h. Each case: 1 CPU/8 GB/3 h. Summary: 1 CPU/4 GB/1 h. Reader:
1 CPU/8 GB/1 h. Every batch is announced before submission. Ordinary jobs only,
no arrays or interactive QOS. Historical per-case placement is retained for
exact comparisons. No running job is cancelled, moved or otherwise disturbed.

The read-only isolated noneditable environment is
`dev/review_depth_balance_20261110/environment`; compute reports Python 3.11.18,
NumPy 2.4.6 and SciPy 1.17.1. Neutral working directories, `-I`, no user site,
PYTHONPATH/PYTHONHOME removal, one BLAS thread and the shared environment read
lock preserve execution isolation. No installations, public/default/test edits,
scientific-criterion changes or frozen-output changes are made.

The independent reader will reconstruct all outside discovery exposures from
source matrices and identity counters, verify exact prior fields and count
conservation, recount all 27 score arrays and nine decisions, and compare 36
deterministically selected cell contrasts against full SciPy NB PMFs integrated
over Gaussian state with explicit tail bounds. Success checks computation only.
No exposure/depth/prior/split/threshold sweep follows a scientific failure.

[Fixed design](../dev/compositional_state_20261111/DESIGN.md),
[bounded model and derivatives](../dev/compositional_state_20261111/model.py),
[shared-noise conditional score](../dev/compositional_state_20261111/joint.py),
[rational/adaptive/gradient checks](../dev/compositional_state_20261111/checks.py),
[workflow](../dev/compositional_state_20261111/study.py),
[executed sources](../dev/compositional_state_20261111/source_manifest.json),
[initial failed source identities](../dev/compositional_state_20261111/attempt1_manifest.json),
[replay correction without tolerance change](../dev/compositional_state_20261111/REPLAY_CONTEXT.md).

[Initial command](../dev/compositional_state_20261111/prepare_submission.json),
[failed equality output](../dev/compositional_state_20261111/job_23207181/output.txt),
[retry command](../dev/compositional_state_20261111/prepare_retry_submission.json),
[independent prerequisite results](../dev/compositional_state_20261111/prepared_v2/NUMERICS.json),
[all input identities/exposures](../dev/compositional_state_20261111/prepared_v2/INPUTS_LOCKED.json),
[replay commands](../dev/compositional_state_20261111/replay_submissions.jsonl),
[complete 27-array qualification](../dev/compositional_state_20261111/qualification_v1/ZERO_SHARE_REPLAY.json),
[case commands](../dev/compositional_state_20261111/case_submissions.jsonl),
[summary command](../dev/compositional_state_20261111/summary_submission.json),
[reader protocol](../dev/compositional_state_20261111/READER_DESIGN.md),
[independent reader](../dev/compositional_state_20261111/reader.py),
[reader source identities](../dev/compositional_state_20261111/reader_manifest.json),
[reader command](../dev/compositional_state_20261111/reader_submission.json).

## Complete scientific result and decision

All nine case fits finish, followed by the complete summary and independent
reader. The fixed screen fails **7/9**. Biological n remains zero. The five
positive-case labels refer to known simulated RNA contrasts; these are not
complete clone-recovery verdicts or independent biological validations.

| Locked case | Predictive advantage | Numerical diagnostic | Smaller discovery mass | Decision |
|---|---:|---:|---:|---|
| original20261031 | -0.776271325 | 0.000740805 | 191.9509 | Fail |
| original20261107 | +6.179664301 | 0.000796963 | 167.8784 | Pass |
| original20261108 | +2.675432046 | 0.000862189 | 183.9726 | Pass |
| original20261109 | +1.860060811 | 0.000523647 | 208.7540 | Pass |
| original20261110 | +1.493709851 | 0.000700798 | 222.9262 | Pass |
| state202611096001 | +5.48e-8 | 0.001742560 | 294.5293 | Unresolved; fails screen |
| state202611096002 | -0.190350421 | 0.000800366 | 176.3481 | Pass |
| state202611096004 | -0.815210867 | 0.001335606 | 129.3809 | Pass |
| state202611096006 | -1.225231278 | 0.001174447 | 292.3576 | Pass |

All component masses exceed 30. The failed positive case is well away from its
numerical boundary. The unresolved state case is neither a measured false-split
event nor a successful resolved negative rejection. It remains failed under the
unchanged declared screen; no rounding, abstention relabeling or threshold change
is used to pass it.

The unbounded shared-noise model scored -0.786648 on 20261031; normalization
changes this to -0.776271, insufficient for positive predictive support. The
20261109 contrast remains positive, changing +1.858881 to +1.860061. No one of
the nine screen decisions changes. This matched experiment does not support
regional mean normalization as a sufficient explanation or repair of the
recovery/state separation failure. It does not establish that compositional
effects are absent generally.

The response law is bounded and its derivatives are checked; that mathematical
property does not supply the missing scientific benefit. Retain the calculation
as a research reference. No w, depth, prior, variance, allocation, basis or
threshold sweep follows. A complete-workflow expansion or new null/power study
for this unchanged failed candidate is **not run**. The original production
2/5 recovery, 14/200 state false splits, detection-power report and 54/1,000
selected-interval coverage are unchanged. No SOTA or genomic membership claim
is supported.

## Independent complete-result checks

Reader23207272 completes 0:0 in 17 seconds. It imports no candidate Observation
implementation. All source/exposure arrays are reconstructed and the declared
score is independently checked using complete SciPy NB PMFs with logistic
fractions, adaptive Gaussian integration and bounded tails.

- All **27** per-cell score arrays sum exactly to the stored values under
  Fraction arithmetic; all order diagnostics and **nine** decisions agree.
- All prior prepared fields are exact, including discovery counts, expectations,
  shapes, depth covariates, q, identifiers and gene masks. D+C matches the
  original source counts in every case.
- The unchanged Philox/domain/binomial construction independently reconstructs
  all outside-region discovery counts. All **4,944** fixed cell observations
  retain their exposure and identity correspondence across the nine cases.
  Fraction reconstruction of the shares differs by at most **1.39e-17**, within
  the declared operation/magnitude rounding diagnostic.
- **36** full-PMF cell contrasts agree with maximum discrepancy **5.85e-14**.
  They use first/last global IDs and min/max depth, without outcome selection.
  Adaptive integration and tail diagnostics remain separate from scientific
  uncertainty and do not certify a global integration enclosure.
- The completed 27-array zero-share prerequisite remains exact under both old
  and new models, with recorded v3/v4 hardware contexts. The initial failure is
  retained and has not been passed by changing a tolerance.

Computational agreement is established for the executed comparisons. Statistical
identification remains conditional on the fitted reference/noise/state law;
biological validation is unestablished. No public package/native/GPU suite is
rerun because this experiment changes only private research code. The prior
validated public source identities remain exact.

## Complete job and preservation record

| Job | Work | Progress | Final state | Elapsed | Process peak RSS, KiB |
|---|---|---|---|---|---:|
| 23207181 | Initial prerequisite | Exact replay failed ██░░░░░░ | Failed 1:0 | 39 s | 397,464 |
| 23207217 | Numeric/data retry | 9/9 inputs ████████ | Complete 0:0 | 14 s | 190,948 |
| 23207237 | Original v4 replay | 24/24 arrays ████████ | Complete 0:0 | 34 s | 283,312 |
| 23207238 | Original v3 replay | 3/3 arrays ████████ | Complete 0:0 | 24 s | 283,776 |
| 23207239 | Complete replay qualification | 27/27 arrays ████████ | Complete 0:0 | 0 s, scheduler | 42,060 |
| 23207253 | original20261031 | 3/3 orders ████████ | Complete 0:0 | 52 s | 260,580 |
| 23207254 | original20261107 | 3/3 orders ████████ | Complete 0:0 | 4m03s | 249,064 |
| 23207255 | original20261108 | 3/3 orders ████████ | Complete 0:0 | 1m47s | 251,516 |
| 23207256 | original20261109 | 3/3 orders ████████ | Complete 0:0 | 1m02s | 244,956 |
| 23207257 | original20261110 | 3/3 orders ████████ | Complete 0:0 | 56 s | 249,904 |
| 23207258 | state202611096001 | 3/3 orders ████████ | Complete 0:0 | 5m40s | 278,540 |
| 23207259 | state202611096002 | 3/3 orders ████████ | Complete 0:0 | 11m44s | 284,012 |
| 23207260 | state202611096004 | 3/3 orders ████████ | Complete 0:0 | 13m18s | 280,964 |
| 23207261 | state202611096006 | 3/3 orders ████████ | Complete 0:0 | 1m28s | 283,416 |
| 23207262 | Complete summary | 9/9 cases ████████ | Complete 0:0 | 5 s | 41,972 |
| 23207272 | Independent reader | 36/36 PMF probes ████████ | Complete 0:0 | 17 s | 205,040 |

Fifteen complete, one failed preserved, zero running or pending in this chain.
Downstream execution is complete; the scientific gate is failed. Every job uses
nslab CPU with the resources listed above. Scheduler sampling reports zeros for
some short jobs; the measured process peaks are reported separately. Maximum
process RSS is 397,464 KiB. No running or unrelated job is disturbed.

All 91 prior validated package/native/data/test identities remain exact. The
284 pre-experiment ledger values and the separate prerequisite checkpoint entry
are preserved. A distinct completed-result entry records the negative result;
it does not rewrite the earlier live checkpoint as if it had already been known.
Source manifests, the failed source snapshot and all earlier frozen results
retain their identities. The full genomic-profile, reliable-clone, spatial and
competitor-superiority goal remains active and unachieved.

[Every case and exact scores](../dev/compositional_state_20261111/summary_v1/REPORT.json),
[independent reader results](../dev/compositional_state_20261111/reader_23207272/results/REPORT.json),
[complete Slurm accounting](../dev/compositional_state_20261111/accounting_complete.tsv),
[process resource records](../dev/compositional_state_20261111/resources_complete.json),
[preservation verification](../dev/compositional_state_20261111/preservation_complete.json).
