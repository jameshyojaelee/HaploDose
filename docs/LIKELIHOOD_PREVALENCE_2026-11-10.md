# Discovery-fitted prevalence loses recovery despite better likelihood

[Observed pattern, not a guarantee] Correctly profiling the discovery mixture
prevalence rejects all3 known state-null splits but also rejects all5 original
recoveries. Every24 historical per-cell predictive array reproduces EXACTLY
on the verified original node. The numerical correction works; the scientific
candidate fails. Reject adoption and full-workflow expansion. Production stays
2/5, all original assertions and5% confidence/power criteria remain unchanged.

## The correction tested and its identifying assumptions

The preceding safeguard fixed pi_s to mean proposal membership q within depth
stratum s. The mixture likelihood identity pi=mean posterior q holds at its own
fitted posterior, not automatically for q from another model. Test whether
re-estimating prevalence under the actual observation model restores information
lost through confirmation. All5 recovery cases and ALL3 known false-split cases
from the latest adjacent-reference candidate are frozen, not a favorable subset.
[Original design](../dev/likelihood_prevalence_20261110/DESIGN.md),
[exact case identities](../dev/likelihood_prevalence_20261110/CASES.json).

Estimand: conditional confirmation RNA predictive advantage of two regional
populations over one continuous state, conditional on the selected proposal.
Both use NB genes, depth-specific location and a Gaussian shared per-cell
regional state. The alternative adds common log contrast delta and stratum
prevalence pi_s. Its two means are alpha_s-log(1-pi_s+pi_s exp(delta)) and that
plus delta. Null, counts, reference noise, offsets, genes, rows, folds, depth
strata, event extent, prior discovery q and support rules remain fixed.

Assumptions include adequate NB means/sizes and latent-state distribution.
Regional transcription, reference mismatch, weak mixture identification and
estimated count-fold dependence remain alternatives. Independent same-cell
DNA would distinguish genotype and validate memberships; independently measured
state/depth and biological replication would test applicability. These8 exposed
mechanism probes have biological n=0 and estimate no biological false-split rate.

Fit ALL nuisance parameters using discovery only at24/48/96 quadrature, write
models BEFORE confirmation counts decode, then evaluate
log p(C|D,mixture)-log p(C|D,continuous), sharing the latent state across folds.
Require exact full=D+C. All5 true cases must have resolved-positive score and
ALL3 known false cases resolved-negative, using the unchanged successive-order
plus rounding diagnostic. This screen is not the5% whole-specimen test. Failure
rejects expansion and threshold tuning; passing would only permit full original
workflow, including merge, followed by fresh complete-workflow error/power.

## Repair numerical failures without relaxing their gates

The first implementation jointly optimized free prevalence logits, arithmetic
means, contrast and variance. Independent full SciPy NB-PMF quadrature agrees
with objective/gradient within4.44e-16 against its1.16747e-9 allowance. An exact
Fraction fixture gives nonzero prevalence score0.04774432884, error6.25e-17.
Nesting, component/cell/stratum symmetry and the variance0 boundary pass.
[Initial numerical evidence](../dev/likelihood_prevalence_20261110/numerics_v1/report.json).

The first8 case attempts stop computationally:4 fail exact old predictions and
4 fail to retain their feasible fixed-prevalence likelihood. These are NOT
interpreted as scientific failures. A separately repeated108 fit on the original
node still fails feasible likelihood, isolating a genuine optimizer issue.

Two controlled old-only replays, without ANY new fitting, isolate the other
failure to node environment:005 reproduces24/24 arrays EXACT;009 reproduces
1/24, with398 changed cell/order entries and maximum7.11e-15. Same wheel,
native hash, package versions and AMD EPYC9254 CPU model; differing microcode
and runtime mappings. This experiment does not identify which particular
runtime/microcode difference caused it. Corrected runs pin the original005;
no tolerance or frozen array changes. This is computational reproduction, not
node selection for favorable scientific outcomes.
[Replay design](../dev/likelihood_prevalence_20261110/REPLAY_DESIGN.md),
[005 control](../dev/likelihood_prevalence_20261110/replay_005_v1/report.json),
[009 control](../dev/likelihood_prevalence_20261110/replay_009_v1/report.json).

Correct the optimizer by profiling pi_s directly. With fixed component
locations/variance, its log likelihood is concave. Scaled component likelihoods
a,b give score sum((b-a)/((1-pi)a+pi b)). Use expm1 for small differences,
compensated sums, boundary score inequalities, then bisection to adjacent
floating-point values or an exact root. No arbitrary root epsilon, prevalence
floor, extra start grid or tuning prior. Flat likelihoods explicitly identify
no prevalence. Use the same7 location/contrast/variance starts and existing
L-BFGS-B settings. Require both profile-anchor inclusion and no loss of its
feasible likelihood. Closed boundaries are the one-component limits of the
unrestricted-logit family. [Optimizer correction](../dev/likelihood_prevalence_20261110/PROFILE_DESIGN.md).

The first profiled numerical fixture incorrectly tested zero nuisance derivative
while holding arithmetic mean fixed; profiling holds component location fixed.
Correct its chain transformation, retaining the SAME allowance and zero-score
requirement. Preserve failed source_v3 and failed numeric output. This is a test
coordinate correction, not a relaxed computational or scientific tolerance.
[Exact correction](../dev/likelihood_prevalence_20261110/NUMERICS_CORRECTION.md).

Corrected numerical checks recover exact prevalence1/6 within2.78e-17, retain
boundary/flat/label behavior and match the independently checked full-PMF
gradient within2.22e-16. Fixed-component nuisance score error5.55e-17.
Five-point gradient errors converge from1.39e-8 to1.50e-12, then3.32e-12 as
rounding increases, below the unchanged1.16747e-9 reference allowance.
[Accepted numerical report](../dev/likelihood_prevalence_20261110/profile_numerics_v2/report.json).

## Complete corrected screen, including every recovery failure

| Fixed case | Cells | Genes | Predictive mixture advantage | Numerical diagnostic | Guard outcome |
|---|---:|---:|---:|---:|---|
|original20261031|508|100|-0.761725259|0.0000568424|reject; loses recovery|
|original20261107|512|93|-2.826582764|0.0000832340|reject; loses recovery|
|original20261108|514|94|-12.821516851|0.0000778681|reject; loses recovery|
|original20261109|502|93|-8.015722997|0.0000537451|reject; loses recovery|
|original20261110|510|94|-21.191226405|0.0000645431|reject; loses recovery|
|state202611097002|600|93|-21.851550928|0.000139066|reject|
|state202611097005|600|93|-19.798655363|0.000105502|reject|
|state202611097007|598|94|-8.067344188|3.02445261|reject|

All8 signs resolve under the prespecified diagnostic. The last false-split
diagnostic is much larger:order sensitivity remains explicit, and this
heuristic is not a formal integration enclosure. The5 decisive wrong-sign
recoveries already fail the screen by margins far beyond their diagnostics.
No extra precision or cutoff tuning could establish a passed screen merely
by resolving another case more tightly.

The corrected alternative fits discovery better than fixed prevalence in all8
cases, while losing confirmation recovery. Several fitted stratum prevalences
hit0/1; original110 ranges from0 to1. [Inference] This is consistent with
overfitting weak composition information or observation-model mismatch, not
proof of a unique cause. Maximizing discovery likelihood does not establish
the genomic explanation or useful power. The all-stratum prevalence parameters
do not create additional independent observations.
[All order-specific results and decisions](../dev/likelihood_prevalence_20261110/profile_summary_v1/report.json).

Every old3 per-cell prediction arrays reproduces EXACTLY in each corrected
case,24 total; counts conserve exactly and no confirmation values fit the new
models. Discovery/joint memberships and effective counts remain in each case's
prediction files; they are diagnostics, not replacement independent confirmation
statistics or cell-level biological truth. No merge change or full-workflow
rescue is claimed. The prior unguarded5/5 research candidate still fails fresh
state-null control. The original3 production failures remain unresolved.

## Every attempt and retained evidence

Research code:
[likelihood_prevalence.py](../benchmarks/hardening/likelihood_prevalence.py),
[profiled_prevalence.py](../benchmarks/hardening/profiled_prevalence.py),
[replay control](../benchmarks/hardening/prevalence_replay_control.py).
This completed result uses immutable source_v4, before the separately designed
pooled-prevalence follow-up. All four source snapshots and failed partials stay.
The working research files may add follow-up modes; reproduction uses source_v4.

| Jobs | Role | Outcome |
|---|---|---|
|23157398|initial independent numerics|complete0:0,17s|
|23157399–23157406|first8 case fits|all failed;4 exact-replay,4 feasible-anchor|
|23157407|dependent first summary|canceled while pending after failed prerequisites|
|23157485–23157486|005/009 old-only controls|both complete0:0,17s each; reported mismatch retained|
|23157506|108 fit on original005|failed feasible-anchor,50s|
|23157573|first profiled numeric fixture|failed nuisance-coordinate assertion,14s|
|23157589|corrected profiled numerics|complete0:0,7s|
|23157605–23157612|corrected8-case screen|all complete0:0,173–359s|
|23157613|dependent complete summary|complete0:0,8s; scientific screen fails|

24 CPU attempts:13complete,10failed,1pending-job cancellation. All are terminal.
No running attempt was canceled; original unprofiled7005 finished with a
feasible-anchor failure after19m02s. All counts above include it. Numerical/
replay/summary jobs1CPU/4GB/1h; case workers1CPU/8GB/3h, innovation/nslab,
CPU/v6&amd, corrected runs pinned005, no arrays/interactiveQOS. Peak process
RSS181348KB, peak scheduler RSS124484K.
[Every exact command](../dev/likelihood_prevalence_20261110/jobs.json),
[cancellation record](../dev/likelihood_prevalence_20261110/cancellations.json),
[accounting](../dev/likelihood_prevalence_20261110/accounting_final.tsv),
[resources](../dev/likelihood_prevalence_20261110/resource_summary_final.json),
[preservation](../dev/likelihood_prevalence_20261110/preservation_final.json).

No production/native/original test/default/dependency/shared-installation or
frozen-result changes. No new full package/GPU suite is claimed:all runs verify
the same isolated wheel previously validated by456 Python/4 Rust tests.
The whole-specimen5% simultaneous confidence/power target and biological
replication requirement remain. The full goal is active and unachieved.
