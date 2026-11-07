# Membership uncertainty and merging

The goal is to distinguish information lost by hard assignments from unsupported
clone distinctions, retaining the original recovery assertions and specimen
error/power criteria. The evidence below uses frozen development simulations;
it adds no biological specimens or independent genomic validation.
Performance findings below are [Observed pattern, not a guarantee].

## What the trace establishes

The three failed cases were replayed from the same count matrices and fitting
metadata, with truth columns excluded from inference. Return-event tracing saved
proposal probabilities, successive profile updates, and merge decisions. Jobs
22849924–22849926 completed without capture errors and reproduced all three
failures. Source, parameters, wheel identity and inputs are recorded under
`dev/membership_20261107`.

For seed 20261108, eight refinement rounds after the accepted chr12 proposal
leave regional profile means 0.03353 and 0.29868: a contrast of about 0.26515
log-rate units. The distinction is then removed by the depth-matched hard-pool
merge test. **The profile has not collapsed before that merge.** This rules out
profile erasure as the immediate explanation for this case. It does not establish
that the fitted memberships or genomic interpretation are correct.

The probabilities used to form event proposals imply substantial contamination
of their hard-assigned pools:

| Seed | Mean high-component probability in low / high assigned pool | Predicted hard-pool log-rate contrast | Proposed latent contrast |
|---|---:|---:|---:|
| 20261108 | 0.251 / 0.777 | 0.1380 | 0.2646 |
| 20261109 | 0.317 / 0.742 | 0.1028 | 0.2443 |
| 20261110 | 0.223 / 0.773 | 0.1393 | 0.2539 |

These are exposure-weighted fitted probabilities, not measured pool purities.
The observation model predicts a pool rate as `(1-q)*theta0 + q*theta1`.
Consequently a hard-pool contrast is a different estimand from the latent
component contrast. The existing HMM's whole-step alternative need not fit the
attenuated pool contrast. Uncertainty in fitted probabilities, nuisance effects
and thinning parameters remains a competing explanation for apparent attenuation.

## Longer refinement is not a recovery correction

The historical refinement loop stops when fewer than 0.5% of MAP labels change.
This is a heuristic, not convergence of probabilities, profiles or nuisance
parameters. An independently solvable contraction fixture demonstrates the
problem: all MAP labels can remain fixed while profile iterates continue moving.

An isolated candidate removed that shortcut and used all eight existing rounds.
Its full-suite job 22850132 recorded **311 passes, one failure and one GPU-only
skip**, plus four Rust passes. The unchanged original thinning assertion fails
because it returns one tumour group. The five-case replay (job 22850142) passes
only **1/5**, retaining seed 20261107 but losing the original seed 20261031; none
of the three failures is recovered. In seed 20261108, normal false calls improve
from 10/90 to 8/90, but recovery still fails. More computation does not resolve
the scientific problem. This candidate and its failed checks are preserved.

The historical default is retained. `SpecimenModel.em` now records the actual
number of refinement rounds and explicitly reports that either MAP stability or
budget exhaustion does **not** certify convergence. Its private `stop_on_map=False`
option supports the controlled budget diagnostic. The contraction test still
requires the exact eight-round result for that mode; no historical scientific
assertion was changed.

Merge logs now report pairwise membership entropy and label the decision as no
resolved difference between assigned pools, **not evidence of genomic equality**.
This records a scientific limitation without changing merge acceptance rules.
The two merges in seed 20261108 have mean pair entropies 0.3779 and 0.4693 nats
(the binary maximum is log(2)). These summarize fitted uncertainty, not measured
misclassification rates.

The retained candidate passes **315 Python tests, one GPU-only skip and four
Rust tests** in compute job 22850847. Three additional targeted statistic tests
pass in job 22851096 (two overlap the full suite). The final unchanged five-case
panel remains **2/5** in job 22850916, with seeds 20261108–20261110 still failing.
The original positive controls are preserved; no new recovery success is claimed.

## Integrating membership uncertainty recovers conditional information

For each locked proposal, the diagnostic compares confirmation counts with two
predictions using identical component rates and nuisance parameters. One assigns
each cell to its MAP component. The other sums the two component likelihoods using
the frozen proposal probabilities, once per cell across the complete region.
Both use the same single-rate NB likelihood, maximized on confirmation data, as
the denominator.

| Seed | Integrated predictive log gain | Hard-assignment predictive log gain | Proposal accuracy | Conditional combined accuracy |
|---|---:|---:|---:|---:|
| 20261108 | 40.993 | -25.024 | 73.66% | 81.19% |
| 20261109 | 7.550 | -29.127 | 69.40% | 74.40% |
| 20261110 | 17.731 | -60.991 | 71.71% | 79.68% |

The last column comes from separately fitting the conditional two-component
model on confirmation counts and updating memberships, then joining truth labels
for scoring among simulated tumour cells in the selected parent. It is **not**
whole-specimen accuracy, an independent membership evaluation, a new pipeline
output, or a full recovery pass. Estimated latent log contrasts are 0.2411, 0.1720
and 0.1872. Their conditional SEs (0.0216, 0.0318, 0.0239) omit uncertainty in
proposal probabilities, nuisance fitting and selection; they carry no calibrated
confidence claim.

The independent full-NB-PMF/Nelder–Mead objective agrees with the gradient fit
within 9.1e-13, and finite-difference checks pass. BFGS nevertheless reports
precision-loss status for two cases; those statuses are retained in `latent_v1`.
The diagnostic evidence supports a better likelihood representation of uncertain
memberships, not a convergence claim obtained by ignoring solver status.

## Complete-search error and power assessment

The experimental statistic is the maximum conditional predictive log gain over
all eligible event proposals encountered in the executed binned workflow,
including proposals later rejected or merged by hard-pool tests. The observer
does not change clone assignments. Calibration reruns initialization, nuisance
estimation, proposal generation, splitting, merging and stopping on each simulated
specimen. Public gene/reference preprocessing and biological samples are outside
this experiment's scope.

A matched control uses hard assignments on **exactly the same proposals**, component
rates, exposure, dispersion and profiled-null denominator. This distinguishes
membership integration from changing which proposals are considered. Raw final
groups remain a comparator. Each statistic has separate fixed NB and
correlated-state calibration banks of 199 draws and 200 evaluations in each of
three conditions. The generator, seeds, rank cutoff of 0.05 and original recovery
assertions are unchanged. Previously examined seeds are paired regression evidence,
not a new independent evaluation.

Before the conditional-search runs, the comparison was locked to five procedures
across three conditions, with one-sided exact simultaneous bounds over all 15
endpoints. Original family-nine summaries are also retained for comparability.
No universal power cutoff is invented. Independent biological replication and
appropriate DNA truth remain necessary for a biological false-split guarantee.

The longer-refinement replay finishes all 998 fits. Its matched NB null has zero
errors, but the state-calibrated gate has 10/200 state-null errors (family-nine
upper bound 10.31%) and 68/200 detections (lower bound 25.71%). The original
procedure had 8/200 and 65/200, respectively. The small detection increase does
not meet the error constraint, and median raw membership ARI remains 0.0764.

The audited even/odd membership-statistic replay gives the following results. Bounds
are simultaneous across the prespecified 15 endpoints; denominators are 200
independent simulation draws per condition, conditional on the realized banks.

| Procedure | Matched NB null: events / upper bound | State null: events / upper bound | Two-clone condition: detections / lower bound |
|---|---|---|---|
| Raw final groups | 0; 2.81% | 184; 96.31% | 200; 97.19% |
| Integrated membership, NB calibration | 0; 2.81% | 47; 32.52% | 49; 16.76% |
| Integrated membership, state calibration | 0; 2.81% | 4; 6.43% | 0; 0% |
| Hard membership, NB calibration | 0; 2.81% | 28; 21.85% | 18; 4.39% |
| Hard membership, state calibration | 0; 2.81% | 5; 7.19% | 0; 0% |

Integration increases both false and true flags with NB calibration. State
calibration suppresses all observed alternative detections for both versions.
Four false events is an observed rate of 2%, but its upper bound does not establish
the requested rate of at most 5%; it does not prove the true rate exceeds 5%.
The original family-nine bound is also above target (6.09%). The negative decision
therefore does not arise solely from adding the matched control to the confidence
family. Neither version qualifies as a recovery correction.

The event statistic is a global RNA flag, not a revised clone caller. It scores
divisive event proposals encountered after initialization; it does not directly
score every initially fitted group contrast. Its zero detection rate rejects this
particular gate as a substitute, **not all possible membership-aware procedures**.
Pipeline memberships are unchanged and remain unsupported. The seed 20261108
normal-cell criterion also still fails (10/90); a diagnostic within a selected
tumour parent cannot repair the normal-versus-tumour assignment.

The first strict audit found different reconstructed latent-mean hashes for the
200 two-clone simulations run on `ne1dc6-010`, while counts, references, dispersion,
truth, genomic resets and state draws were identical. The same software versions
were recorded. The failed audit and first replay remain preserved. Only that
shard was replayed on the original `ne1dc4-063` host, with a fresh output and
all identity checks intact. Audit job 22852111 verifies all 998 input records in
each experiment, independently recounts decisions and inverts binomial tail
probabilities for the bounds. The membership observer leaves all 998 raw
assignments and posterior arrays exactly unchanged from the original procedure.
The rejected longer-refinement candidate changes 166 of those assignment sets.

The original-host replay changes proposal eligibility in two alternative
specimens (replicates 150 and 167): two integrated-membership decisions swap,
leaving their total at 49, and one hard-control acceptance disappears (19 to 18).
These are three CPU-sensitive experimental decisions, not merely different
printed decimals. The immediate discrepancy is presence versus absence of an
eligible event proposal; its upstream numerical cause is not localized here.
Both versions remain preserved. Cross-hardware numerical certification is not
established. The negative error/power conclusion is unchanged.

The original benchmark uses even/odd bins, whereas the three recovery failures
use count thinning. A separate, mode-matched comparison therefore changes only
`split_mode` to `thin`, keeping the same counts, seeds and decision criteria.
Its design additionally declares simultaneous bounds across all 30 endpoints
from both modes. All 998 thinning fits and the independent mode audit completed.
The key state-calibrated results under that joint confidence family are:

| Mode and membership treatment | State-null flags / 200 | Joint upper bound | Alternative detections / 200 |
|---|---:|---:|---:|
| Even/odd, integrated | 4 | 6.88% | 0 |
| Even/odd, hard | 5 | 7.66% | 0 |
| Thinning, integrated | 10 | 11.26% | 0 |
| Thinning, hard | 7 | 9.15% | 0 |

In thinning mode, NB calibration gives 36 false flags and 53 detections with
integration, versus 16 false flags and 19 detections with hard assignments.
Increased detection accompanies increased false distinctions; state calibration
again removes all observed detection. Every method has zero matched-NB-null
events, which does not establish robustness to the correlated-state condition.
The complete endpoint table and paired comparisons are in
`dev/membership_20261107/mode_evidence_v1/report.json`. Paired hard-versus-integrated
tests use exact binomial tests on discordant specimens, with Holm correction
over 12 comparisons; they do not replace the error-bound and power criteria.

Job 22852791 checks identical inputs for all 998 thinning specimens, independently
recounts ranks and decisions, and inverts binomial tails for the bounds. Three
prespecified controls (evaluation replicate zero in each condition) rerun the
thinning fit without the observer and reproduce posterior arrays exactly.
No inference-rule replacement is justified by these experiments. The three
complete recovery failures remain unresolved.

## Identification and decision limits

Computational agreement establishes the implementation of a likelihood or count;
it does not identify genomic cause. A continuous regional expression state can
produce the same RNA rate distribution as the proposed CNA components. A global
rank rejection distinguishes the observed RNA from a specified calibration law;
it neither certifies every proposed pair nor validates cell memberships. Both
membership-statistic procedures therefore retain unsupported membership status.

The proposed correction is falsified as a general clone solution if it exceeds
the prespecified false-split bound, merely suppresses both true and false
contrasts, or gains apparent power only from unmatched proposal searches. A
positive conditional likelihood on the three inspected cases is insufficient to
promote it. The original failed outputs and all negative experiments remain in
the evidence record.

## Preservation

Compute job 22851623 verifies 340 initial source files, 340 files in the previous
final source, nine recorded historical-result hashes, four original
criteria/generator files, all 20 recorded inputs from the five-case panel, and
four wheel/native pairs. All 22 current package source/data files match the
tested candidate. This bounded audit does not claim an exhaustive inventory of
historical raw data. Its record is `dev/membership_20261107/preservation_audit.json`.
