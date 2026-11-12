# Fixed working RNA dosage against continuous-state prediction

The fixed working dosage model fails the unchanged screen. Four of five
original cases favor the discrete mixture. All four known state false
splits favor the continuous-state model. The required 20261109 contrast
improves from −3.396705 to −1.955938 log-predictive units, with the sign
still resolved against the mixture. All 27 historical per-cell arrays
and corresponding scores reproduce exactly. The candidate stays in
research and the production recovery count remains 2/5.

The depth-conditioned free-contrast mixture recovered the expected
predictive sign in four of five original cases and all four exposed
state false splits. Original 20261109 still favored continuous state.
This experiment tests a specific source of predictive loss: estimating
the separation of two mixture components from the same weak proposal
reads that select the candidate. It fixes the magnitude to the existing
working dose step and preserves the continuous-state alternative.
[Prior screen](DEPTH_STATE_PREDICTIVE_2026-11-10.md),
[locked design](../dev/fixed_dosage_state_20261110/DESIGN.md).

## Scientific and numerical contract

The estimand is conditional RNA log-predictive gain of the selected
two-component model over a depth-conditioned Gaussian regional state.
Both use the same shared latent state across discovery and confirmation
and allow continuous within-component variance. The fixed alternative
uses the archived native one-copy step, `beta*log(1.5)`. Its sign is
selected on discovery reads from the two possible orientations. All
other fitted locations, variance starts, depth-stratum prevalences,
cells and gene identities follow the prespecified design. The null
and the previous free-contrast model reuse their exact saved discovery
fits. A separate arm with the same native noise model checks that
the original per-cell scores and aggregate decisions reproduce.

The largest unresolved biological alternative is regional transcription
producing the same observations as a genomic dosage change. Independent
corresponding DNA at the claimed resolution can distinguish these cases.
The count folds also rely on estimated noise; the computation does not
establish their sampling independence. The sign gate must retain all
five original cases and reject all four exposed state false splits. A
pass would permit the original full-workflow recovery and false-split /
power gates. These selected nine cases supply no biological replication.

Numerical job 23188145 passed the inherited independent PMF and
quadrature checks and the new reduced-parameter checks. Direct exact
rational arithmetic verifies the zero-variance objective and gradient
to 2.50×10⁻¹⁶ and 5.55×10⁻¹⁷, within an operation-derived 2.27×10⁻¹³
budget. Independent full-PMF adaptive integration with centered finite
differences verifies the positive-variance gradient to 5.11×10⁻¹³.
The difference across successive extrapolations plus rounding gives
a 4.67×10⁻¹⁰ numerical diagnostic. Component exchange is preserved.
[Numerical checks](../dev/fixed_dosage_state_20261110/numerics_v1/report.json).

## All nine fixed cases

| Case | Cells | Free contrast | Fixed contrast | Fixed numerical diagnostic | Required sign passes |
|---|---:|---:|---:|---:|---|
| Original 20261031 | 508 | +4.348085 | +2.786425 | 0.000052 | Yes |
| Original 20261107 | 512 | +4.713200 | +3.770467 | 0.000066 | Yes |
| Original 20261108 | 514 | +4.527263 | +2.553723 | 0.000066 | Yes |
| Original 20261109 | 502 | −3.396705 | −1.955938 | 0.000037 | No |
| Original 20261110 | 510 | +2.863255 | +5.808933 | 0.000054 | Yes |
| State 202611096001 | 599 | −2.277543 | −4.499497 | 0.000095 | Yes |
| State 202611096002 | 600 | −1.836996 | −3.989898 | 0.000077 | Yes |
| State 202611096004 | 599 | −4.167723 | −4.564059 | 0.000090 | Yes |
| State 202611096006 | 600 | −10.889296 | −11.551857 | 0.000069 | Yes |

Every fixed sign resolves against the original successive-quadrature and
rounding diagnostic. The fixed alternative makes all four state scores more
negative. It also makes three of four passing original scores smaller. The
required 20261109 remains incompatible with this predictive rule. The
fixed contrast cannot serve as a sufficient correction for membership loss
or establish the 5% specimen guarantee. The full workflow and larger
error/power runs were not run for this candidate after the screen failed.
[All scores and decisions](../dev/fixed_dosage_state_20261110/summary_v2/report.json),
[summary execution manifest](../dev/fixed_dosage_state_20261110/summary_v2/run_manifest.json),
[20261109 locked discovery fits](../dev/fixed_dosage_state_20261110/case_03_v1/DISCOVERY_MODELS_LOCKED.json).

## Execution

All jobs use nslab CPU, one CPU each and a one-hour limit. The numerical
prerequisite and the summary use 4 GB each. The nine case jobs request
8 GB each, depend on successful numerics, and read the existing immutable
case and source snapshots. The final summary depends on all nine cases.
The validated read-only installed wheel matches the original model
comparison. No shared installation, frozen output or threshold is edited.
Eight first-pass cases and the numerical prerequisite finish successfully.
Case 06 on node003 stops at its exact old-array comparison, with 22/600
differences at the first quadrature order and maximum absolute difference
4.44×10⁻¹⁶. The dependent summary enters a terminal canceled state after
the dependency fails. The separate cancellation command reports an already
invalid job ID. Case 06 is retried with the original native-runtime setup.
The pending retry is moved from occupied node009 to node010. It completes
with all three original arrays and scores exact. The exact check is
unchanged. The replacement summary uses the same scorer and an explicit
inventory of the nine finished cases. Eight cases take 6–83 seconds each,
the retry takes 49 seconds, and the summary takes 4 seconds. The maximum
measured process peak is 474,412 KiB. Across 13 submitted jobs, 11 finish
successfully, one fails the parity check and one dependent job is canceled.
[Case identities](../dev/fixed_dosage_state_20261110/CASES.json),
[executed source](../dev/fixed_dosage_state_20261110/study.py),
[initial submissions](../dev/fixed_dosage_state_20261110/jobs.jsonl),
[failed parity check](../dev/fixed_dosage_state_20261110/analysis_23188152/output.txt),
[retry command and cancellation response](../dev/fixed_dosage_state_20261110/retry.jsonl),
[pending placement update](../dev/fixed_dosage_state_20261110/retry_placement.json),
[replacement summary command](../dev/fixed_dosage_state_20261110/summary_retry.json),
[all Slurm accounting](../dev/fixed_dosage_state_20261110/accounting_complete.tsv).

The later discovery-prevalence replay qualifies the numerical environment more
precisely. The eight first-pass cases use copied runtime003, while the successful
case06 retry uses the native node010 interpreter/library launch. Moving the
copied launcher to node010 does not reproduce case06. Conversely, native node010
does not reproduce case00. Matching each case's archived qualified runtime
preserves all 27 arrays and original scores exactly; all new discovery inputs
are byte-identical across these launches. The earlier node association alone
does not isolate this library difference.
[Bounded diagnosis and unchanged exact gate](PROFILE_MEMBERSHIP_PRIOR_2026-11-10.md).
