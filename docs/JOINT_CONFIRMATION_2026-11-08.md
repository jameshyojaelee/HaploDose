# Shared-membership confirmation and matched bulk operators

The goal remains accurate regional RNA dosage and appropriately supported CNA
and clone inference, with independent evidence at the claimed resolution.
This checkpoint improves an experimental recovery path from 3/5 to 4/5 without
changing any original criterion. It does not complete the scientific plan or
establish competitor superiority, biological clone validity, or error control.

## Recovery and the uncertainty calculation

### Unordered component labels

Job22976607 tests whether arbitrary component names hide a distinction between
two RNA-rate trajectories. It uses the same three saved events and full-count
model as the failed ordered joint probe, without oracle memberships. The
estimand is now an **unordered pair of trajectories**. Component1 is defined as
the higher-rate component at the first nonneutral bin; one exchange applies to
the entire trajectory and all cell labels. Opposite signs at later loci survive.
The all-neutral trajectory has unidentified labels and contributes0.5 to the
arbitrary membership probability. Anchored cell-specific priors are rejected.

| Seed | Canonical core | Minimum decoded-state posterior | All-neutral trajectory probability |
|---|---|---:|---:|
| 20261108 | bins209–217 | 0.995820 | 2.01e-24 |
| 20261109 | bins209–217 | 0.955557 | 3.68e-6 |
| 20261110 | bins209–217 | 0.985139 | 8.98e-17 |

These pass the existing0.95/five-bin rule for the new estimand. They do not turn
the original ordered tests into passes. Evidence and neutral mass are unchanged;
the label operation is applied to whole trajectories before marginalization and
MAP selection, not by taking independent absolute states at each bin.

Independent SciPy full-NB enumeration over states **and** cell assignments agrees
within3.34e-16 for state probabilities and2.23e-16 for memberships on three small
fixtures. Cell permutation, component exchange, opposite-sign trajectories and
refusal of anchored priors pass. Numerical allowances derive from operation
counts and likelihood scales; no decision lies within its rounding allowance.
The audit verifies425 frozen baseline files,66 unchanged core/test/native files
and26 inputs. Runtime23.93s, peak RSS193696KiB, CPU backend, corrected isolated
wheel. Evidence: `dev/label_invariance_20261108/probe_v1`.

[Inference] Unidentified component names explain the failed directional support
in this joint diagnostic. Whether usable memberships survive profile refitting
and merging remains untested here. This is selected synthetic-event evidence,
not independent confirmation, complete-workflow recovery, error calibration or
biological validation. No package/default or original scientific test changes.

The new research model retains each cell's NB likelihood and a single latent
membership shared across all bins in a candidate interval. It integrates that
membership conditional on every dosage-state trajectory. The state grid, prior,
transition penalty, boundary transitions and every-bin posterior rule are
retained. Thinning uses the actual adjusted minimum of **five bins** and
posterior **0.95**. Some earlier task STATE notes said three bins; executable
thinning runs used five. The original `ProfileParams` minimum of three is not
the thinning-adjusted minimum.

Exact summation is available for short intervals only, with a declared cap of
one million state paths. Longer proposals explicitly retain the previous
confirmation procedure. This is a research hybrid, not a scalable replacement
or a new inference default.

Independent enumeration over both memberships and states, using SciPy NB PMFs,
agrees with the optimized marginal calculation. Boundary, label permutation,
cell permutation, pure-membership, subset component-rate preservation and
variational lower-bound checks pass on compute nodes. Numerical allowances are
tied to float64 accumulation scales and independent reference discrepancies.

The mean-field approximation is **not accepted for decisions**. On seed
20261109, exact nonneutral probabilities range from about 0.54 to 0.73, while
the converged approximation gives approximately 0.999 throughout the region.
Tighter convergence changes probabilities by only 5.9e-11. This is approximation
error, not an optimizer tolerance problem. Exact inference rejects the region
at the unchanged 0.95 rule. Seeds 20261108 and 20261110 retain eight-bin exact
cores; their uncertain ninth bins are excluded.

| Experimental path | Original assertions passed | Interpretation |
|---|---:|---|
| Frozen merge-weighted control | 3/5 | Existing reference, no new default |
| Exact initial confirmation only | 3/5 | Final posteriors identical to control on all five cases |
| Exact initial and depth-matched confirmation | 4/5 | Recovers seed 20261110; other four outputs identical to control |

The recovered case has assignment accuracy **0.812**, normal false calls
**6/90**, tumour recall **500/510**, exactly two tumour groups and the required
chr12 exported difference. All original assertions are executed without edits.
In the first proposal its exact log evidence is 12.192 before depth matching
and 6.910 after matching; the pooled procedure retains no segment at either
stage. Replacing only the first stage therefore does not recover the case.
The existing weighted merge retains the recovered distinction through export.

Depth matching uses the original eight strata, equal integer counts from each
proposal group, and a fixed hash ordering of cell identities. Selection does
not use confirmation counts, probability magnitude or truth. Subsetting retains
component rates normalized on the full proposal population. This avoids treating
fractional cell weights as an exact NB sampling likelihood. The original minimum
cell support, gene veto, ranking statistic, initialization and merge procedure
remain active. The ranking statistic is not newly claimed to be calibrated.

Seed 20261109 still fails the initial exact confirmation. A separate postfit
oracle diagnostic restores its nine-bin core when supplied true memberships,
with fitted exposures and the original state grid retained. Hard proposal labels
fail all three diagnostic cases. [Inference] Membership information is limiting
under this conditional model; this does not prove that the observed RNA is
intrinsically nonidentifiable or that an observable correction exists. Model
centering follows the supplied memberships, so the oracle intervention includes
that specified coupling. Oracle labels never enter pipeline or calibration fits.

## Separating membership information from rate centering

Job 22972350 separates the two changes coupled in the original oracle: cell
membership probabilities in the marginal likelihood and exposure-weighted
component-rate centering. All four combinations run on each of the same three
saved events. Counts, exposures, dispersion, state grid, prior and the original
0.95 posterior/five-bin rule stay fixed. Off-diagonal combinations are controlled
conditional-model interventions, not newly fitted self-consistent mixtures.

For seed 20261109, original soft memberships yield no qualifying core with either
centering: minimum nonneutral probabilities are 0.5424 with original centering
and 0.1840 with oracle centering. True memberships with **original** centering
restore all nine bins, with minimum nonneutral probability 0.9983. Changing both
gives 0.9997. The two other events retain their original eight-bin cores with
either centering and soft memberships; true memberships restore nine-bin cores.

[Inference] Membership information alone is sufficient for the conditional oracle
rescue; changing centering alone is not. This resolves the earlier attribution
ambiguity without identifying an observable correction. It does not show that
RNA contains enough information to recover those memberships, validate posterior
coverage, or justify feeding confirmation-derived memberships back as independent
evidence. Oracle inputs remain outside inference and calibration jobs.

Six diagonal results reproduce the original posteriors, log evidence and cores
exactly. Independent full SciPy NB sums over all states and memberships pass all
16 small-fixture factor/boundary combinations and cell-permutation checks. Maximum
evidence and posterior discrepancies are 1.12e-15 and 6.67e-16. The job completes
in 43 seconds, with 132,956 KiB scheduler peak RSS; 422 baseline files, 66 unchanged
package/test/native files and 21 input hashes are verified. All twelve outcomes,
the prespecified interpretation and execution identity are retained under
`dev/membership_centering_20261108`. No package correction or criterion change
follows from this nondeployable oracle.

## Independent proposal optimizer follow-up

Job22967413 tests whether a poorer local mixture solution explains the remaining
membership loss. All three original saved events, including the two cases now
recoverable by the experimental path, receive ten fixed starts and two solvers:
L-BFGS-B with independently derived gradients and Nelder–Mead with full SciPy NB
probabilities. Proposal counts alone select solutions; confirmation and simulated
labels do not select starts or fits. Mixing-weight bounds and the observation
model remain unchanged.

All60 fits converge. The greatest apparent improvement over the saved solution
is3.65e−12 log-likelihood units, below both independent scalar-reference discrepancy
(maximum1.27e−11) and the derived full-PMF accumulation budget. Saved solutions
are retained in all three cases. Package null-relative objectives agree with the
independent PMFs within1.39e−12; analytic gradients agree within1.54e−15 per cell.
Finite-difference step sweeps show truncation-error decline followed by rounding
error, with best discrepancies below1e−9 per cell. These are numerical checks,
not scientific decision thresholds.

Exact confirmation remains identical: seeds20261108 and20261110 retain their
eight-bin cores; seed20261109 has no qualifying core. The search gives no
justification for an optimizer correction. It is a finite multistart investigation,
not a global-optimum certificate or proof of biological nonidentifiability.
The known membership and observation-model limitations remain.

The job completes in45s on an Intel v4 CPU node, with4 CPUs and8GB requested.
It verifies406 baseline files and65 unchanged package/test/native files.
Sources, all solver attempts, scalar checks and exact posteriors are preserved
under `dev/membership_optimizer_20261108/probe_v1`. No package default or criterion
changes. The investigation does not add biological replication or an error bound.

## Joint-count information probe: no directional recovery repair

Jobs 22974351 and 22974536 test a different conditional calculation on the same
three saved events. Proposal and test counts enter jointly once, with doubled
exposure and NB size. Each cell starts with the same proposal-fitted population
mixing probability; proposal-derived cell memberships are not reused as a prior
for those same counts. Fitted nuisance and the original state grid/prior stay
fixed. This is a selected-event diagnostic, not independent confirmation or a
complete recovery workflow.

| Seed | Minimum nonneutral probability | Minimum decoded directional probability | Original qualifying core |
|---|---:|---:|---|
| 20261108 | 0.9958 | 0.6896 | None |
| 20261109 | 0.9556 | 0.6255 | None |
| 20261110 | 0.9851 | 0.9633 | Nine bins |

The unchanged rule requires the **decoded directional state** to have posterior
at least 0.95 at every retained bin, with at least five contiguous bins. Adding
the probabilities of both nonneutral directions would answer a different
question and does not satisfy that rule. [Inference] Removing cell-specific
membership information leaves substantial direction/label ambiguity in the first
two conditional mixtures. This does not establish that the full observations lack
information about unordered groups. A future joint search would have to define
that label-invariant estimand explicitly and establish its whole-workflow error
and power, rather than relabel unsigned support as directional confidence.

The proposed joint calculation does not rescue the remaining event and loses
the original eight-bin directional core in seed20261108. Per the frozen design,
do not launch a large candidate sweep from this probe. The original conditional
procedure and every historical failure remain intact; no default is promoted.

Numerically, the full-count component log likelihood agrees with the sum of the
half likelihoods to 3.56e-15 before integrating one shared membership and state
trajectory. Exact posteriors agree to 4.89e-15. These identities condition on the
NB working law, not its adequacy for aggregated genes or estimated thinning
nuisance. Sixteen independent full-PMF fixtures pass and all three original
conditional results reproduce exactly.

The first probe used nonneutral probability in its numerical-boundary warning,
although its actual core decisions already used the correct directional rule.
Version2 corrects that diagnostic flag, exports the decoded path and directional
support, and reproduces all12 original posteriors, evidence values and core
decisions exactly. Both versions are preserved; version2 completes in33s with
134,668KiB peak RSS. It verifies424 baseline files,66 unchanged package/test/native
files and22 inputs. No decision lies within the derived rounding allowance of
the actual directional threshold. Evidence is in
`dev/joint_count_information_20261108/probe_v2`.

## False-split and power challenge: completed

The earlier binned experiment has 60-bin chromosomes and 30-bin half-chromosome
proposals. Neither enters the short-event exact path. Repeating it unchanged
would not test the new mechanism. Its failed error/power results remain intact.

A separately locked resolution challenge retains the original generator,
600 cells, 240 bins, counts law, state and CNA effects, seed keys, 199 calibration
specimens per null and 200 evaluation specimens per condition. Chromosome resets
are every 20 bins, producing ten-bin half-chromosome proposals. The old 60-bin
state/CNA region spans three such blocks. This is an artificial resolution
challenge, not a physical chromosome model or new biological replication.

The merge control and exact-matched candidate each run 998 complete binned
workflows. Initialization, nuisance fitting, proposal selection, confirmation,
merging, stopping and final-group scoring run identically in calibration and
evaluation. Raw groups, accepted distinctions, membership ARI and exact-two
counts are retained. There are no gene or allele observations in this generator;
it does not validate those channels or the public preprocessing pipeline.

The 18 new endpoints join the 27 previous endpoints in a **45-member** family
of simultaneous one-sided exact bounds, with tail alpha 0.05/45. The specimen
false-split target remains 5%; no count, alpha or seed adjustment is permitted
after outcomes. Earlier bounds are recalculated only in the new report, alongside
unchanged historical reports. Power is assessed simultaneously; a global
distinction does not certify its memberships. Bounds condition on the realized
calibration banks and cannot replace independent biological specimens.

Jobs 22962400–22962409 and summary22962535 completed. Both experimental paths
have the same marginal error and detection counts:

| Condition and decision | Events / specimens | Simultaneous bound |
|---|---:|---:|
| Matched NB null, all three rules | 0/200 | Upper 3.34% |
| Correlated-state null, raw or NB-calibrated | 187/200 | Upper 97.62% |
| Correlated-state null, state-calibrated | 12/200 | Upper 12.93% |
| Two CNA populations, raw or NB-calibrated | 200/200 | Lower 96.66% |
| Two CNA populations, state-calibrated | 88/200 | Lower 33.36% |

The exact candidate loses one and gains one state-calibrated detection; all
other decisions agree. Paired Holm p-values are 1. Exactly two raw tumour groups
occur in only 1/200 alternatives in either path; median raw ARI is 0.0750 for
control and 0.0752 for the candidate. Detection is not reliable membership
recovery. The 88 detections cannot be attributed to exact confirmation relative
to the previous 65: the genomic reset changed and the new control also has 88.

The exact path is active: its five shards execute 104, 668, 115, 649 and 826
exact likelihood calls. Thus absence of an aggregate benefit is not explained
by a completely inactive intervention. **Do not promote this hybrid.** Its
4/5 recovery result does not overcome failure of the declared error/power test.
Failure of an upper-bound criterion alone does not prove a true rate above 5%.

Audit22962536 verifies all 1,996 generated inputs, 3,600 decisions and 45 bounds
using independent binomial-tail inversion. Its initial coverage omitted direct
summary-count reconciliation. Audit22963491 adds independent recounting of all
45 endpoints and nine paired discordances and passes. Both audits are retained
under `null_audit_v1` and `null_audit_v2`; the simulations were not rerun.

A secondary postfit diagnostic (job22963381) distinguishes intrinsic overlap
from algorithmic membership loss. For the generator's Gaussian state amplitude,
known effect and equal independent class prior, observing the noiseless amplitude
gives expected Bayes accuracy Phi(delta/(2 sigma)) = 0.686714. Independent numerical
integration agrees. The 200 saved alternatives give 0.684892 for that oracle
classifier; a separate known-balanced-count oracle gives 0.685033. These are
model-conditional expected accuracies, not bounds on every realized specimen or
claims of biological identifiability. Cell identities carry no oracle information.

One-to-one group matching, penalizing extra groups and unassigned cells, gives
mean accuracy 0.438292 for control and 0.439617 for the candidate. The candidate
assigns 99.61% of cells; applying the specimen gate as abstention reduces mean
coverage to 43.81% and correct assignments across all cells to 19.10%. The gap
to the noiseless oracle warrants tracing initialization and merging. Neither
oracle labels nor amplitudes entered inference, calibration or parameter choice.

## Bulk accuracy under a shared operator

The raw COAD/READ development panel reuses nine frozen fits: 20 reference donors,
24 calibration patients and 24 previously inspected reserved patients. It is
not a new independent evaluation. Reconstructed sample order, DNA values,
loss/neutral masks and common callable support match the frozen panel exactly.
All original patient-score fields reproduce to text precision except the
separately investigated float32 Pearson reductions.

Naive log2 profiles are converted to natural-log units before sharing a penalty
with HaploDose ratios. The exact least-squares partition operator chooses each
method's own boundaries using RNA only, the same naive first-difference noise
scale, penalty, arm/gap breaks and minimum three approximately 10-Mb windows.
This coarse operator is a diagnostic, not a proposed default. DNA enters scoring
and calibration thresholds, never boundary selection.

The original five methods and six contrasts are retained. Two partitioned
methods and two contrasts yield 24 endpoint comparisons. The analysis uses
10,000 nested resamples of the eight existing donor-reference refits,
calibration patients and reserved patients. Holm tests and Bonferroni basic
bootstrap intervals use the full family. Eight reference draws do not establish
coverage of these approximate intervals.

| HD partition minus naive partition | Paired mean effect | Simultaneous interval | Evaluable patients |
|---|---:|---:|---:|
| Loss recall | −0.0611 | [−0.5129, 0.2584] | 13/24 |
| Spearman correlation | −0.0291 | [−0.0556, 0.0132] | 23/24 |
| Pearson correlation | −0.0588 | [−0.1016, 0.0178] | 23/24 |

No shared-operator advantage is established. The naive partition's mean
Spearman is 0.6713, versus 0.6422 for the HD partition and 0.6699 for existing
HD segmented output. Those marginal means do not establish equivalence or
superiority of the strongest baseline against every package output. The
specific HD-segmented versus naive-partition contrast was not in this locked
primary family and remains an additional descriptive/replication question.

Mean loss-side FPR is 0.0412 for the naive partition and 0.0423 for the HD
partition. Their recalls are 0.3997 and 0.3386 respectively. All seven methods
have complete classified-truth coverage on common support, but endpoint
eligibility differs: 13 patients have loss truth, 23 have evaluable correlations,
and precision is undefined without calls. These denominators are retained.
There is no justified default promotion from the operator experiment.

The first bulk attempt, job22962209, stopped on the Pearson reproduction check.
The historical scorer uses float32 window arrays; float64 text-roundtrip
tolerance was inappropriate for its reduction arithmetic. Each old and rebuilt
correlation was checked against an independently derived forward-error enclosure
covering mean, product, sum, norm and division rounding. New correlations use
extended-precision centered products, independently checked with binary64
`math.fsum`. Maximum old/rebuilt difference was 1.1921e-7; independent refined
calculations differ by at most 2.2205e-16. The failed run and original scores
remain frozen. No scientific threshold was changed to accommodate this failure.

## Evidence, implementation and limits

All new implementation is under `benchmarks/hardening`; production package,
native code and original tests are unchanged. The numerical model and stage
hooks are in `shared_membership_hmm.py` and `shared_membership_stages.py`.
`shared_membership_probe.py`, `shared_membership_panel.py` and
`shared_membership_oracle.py` supply numerical and recovery evidence.
`raw_operator_panel.py` and `pearson_reference.py` implement the bulk diagnostic.
`exact_stage_null.py` runs the activated error/power challenge.
`joint_confirmation_audit.py` and `exact_null_audit.py` perform separate audits.

Evidence is under `dev/joint_confirmation_20261108`, including the locked
protocols, initial dirty baseline, executed snapshots, isolated-wheel identity,
commands, source/input hashes, per-job resources, failure logs and output-to-score
links. Completed jobs are 22961863 (numerical probe), 22962126/22962127 (recovery),
22962278 (bulk), 22962460 (oracle) and 22962501 (audit). The recovery candidate
took 6m46s with process peak RSS about 1.38GiB; bulk took 6m34s with about 621MiB.
The full unchanged package suite was not rerun; its last recorded validation is
319 Python passes, one GPU-only skip and four Rust passes in job22855685.

Audit22962501 verifies 367 baseline and 367 prior-final files, 65 unchanged
package/test/native files, nine recorded historical results, 20 panel inputs,
eight build artifacts and all 24 paired bulk effects. This is a bounded historical
inventory, not a claim to have hashed every historical dataset. External bulk
truth files have a supplemental post-execution hash inventory; reconstructed
DNA/support and old scores were checked before new comparisons.

Computational validation supports this exact conditional model and these
recorded executions. Statistical identifiability remains conditional on fitted
memberships, reference/nuisance information and the observation model. Biological
validation is not established. Routine interval coverage, independent clone error
and power, scalable long-event inference, corresponding DNA membership truth,
spatial validity and a sufficiently replicated unopened biological panel remain
open. More simulated specimens cannot close those biological evidence gaps.
