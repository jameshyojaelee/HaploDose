# Ordered RNA amplitudes under a common count objective

The known-law count oracle passes the original32 RNA and64 both-target
diagnostic comparisons, while the fitted five-structure selector still
fails low-depth minority comparisons. We therefore test a parsimonious
ordered RNA contrast under the same count likelihood. No true profiles,
memberships, population fraction or latent state draws enter fitting.

## Scientific contract and competing explanation

The estimands remain posterior expected regional RNA log odds and the
working response with its continuous expression state set to zero.
Observations are integer multinomial counts conditional on library total,
with an optional independent standard-Gaussian state and one/two RNA
populations. The new contrast has three fitted RNA levels: zero, one
negative amplitude, and one positive amplitude. These are not DNA copy
labels. Shared amplitudes, ordered regions and state independence are
identifying restrictions. Regional transcriptional programmes and
heterogeneous dosage response remain strong competing explanations.
Only corresponding independent DNA at the claimed resolution can support
genomic attribution. Biological replication n=0 here.

The count model and whole-cell heldout K selection are unchanged. A proper
integrated-break path prior replaces arbitrary pair-grid constraints. Its
break probability is integrated under Uniform[0,1], with uniform initial
state/alternative transitions. This is a declared working regularizer,
not a learned genomic transition law or calibrated uncertainty model.
It introduces no selected penalty or copy-state grid. Supplied genomic
discontinuities reset the path prior.

The continuous step optimizes the count likelihood with the path fixed.
An exact dynamic program maximizes its old-posterior multinomial minorizer
plus the same path prior with amplitudes fixed. Jensen and log(x)<=x-1
establish the bound. This does not transplant a Poisson sampling model.
Every update must improve the common count-MAP target within an explicit
forward arithmetic budget. Global nonconvex optimality is not certified.
Fixed native PCA/two-means starts and exact nested parents are retained;
there is no new start/rank/amplitude/threshold search.

Here,count-MAP means unpenalized continuous-parameter likelihood optimization
and a working log-prior penalty on the discrete path. It is not a normalized
joint Bayesian posterior over amplitudes,baselines and paths. Reported latent
expectations condition on fitted parameters;they do not integrate parameter
uncertainty,define genomic membership probabilities or provide validated
confidence intervals. Retaining a feasible parent preserves its valid
prediction and objective value,not a new local/global optimality certificate.

## Mathematical gate completed

Job23222260 completes0:0 on nslab CPU. Exact rational prior normalization
and exhaustive DP comparisons pass for lengths1–5 and the(2,3)block reset,
with zero and seeded random emissions. Independent complete multinomial
probabilities verify minorizer coefficients, contact, raw-gauge identity,
54 alternative likelihood bounds, and the path update. Independent
complex-step and halved-step projected gradients agree, maximum gap
2.776e-16. Population-label/state-sign invariants pass. Allowances follow
operation counts and observed magnitude; no scientific cutoff changes.
This establishes these computational identities, not profile accuracy.

## Mechanism gate

The original eight fixed mechanism cases are next, with all original
heldout selection, numerical diagnostics and mechanism criteria preserved.
The reader independently reconstructs full PMFs, points, actual counts,
all intermediate count/MAP gains and fitted path parameters. All eight
prediction locks must complete before any target opens. Passing this
gate authorizes the original290-case32/64 regression; it does not justify
promotion. A failed mechanism gate blocks expansion. Full-workflow false
splits, detection power, uncertainty and biological validation remain
separate gates. Public methods,30-cell support, frozen outcomes and prior
ledger values remain unchanged. Results pending.

[Fixed design](../dev/ordered_amplitude_joint_20261111/DESIGN.md),
[count minorizer and DP](../dev/ordered_amplitude_joint_20261111/ordered.py),
[common-objective fitter](../dev/ordered_amplitude_joint_20261111/fitter.py),
[independent qualification](../dev/ordered_amplitude_joint_20261111/qualify.py),
[mathematical result](../dev/ordered_amplitude_joint_20261111/qualify_23222260/results/REPORT.json),
[mathematical snapshot](../dev/ordered_amplitude_joint_20261111/SNAPSHOT.json),
[mechanism reader](../dev/ordered_amplitude_joint_20261111/reader.py).
[Completed independent eight-case reader](../dev/ordered_amplitude_joint_20261111/reader_23222317/results/REPORT.json).

## Mechanism result:pass with low-depth failures retained

All eight writers23222309–23222316 and reader23222317 complete0:0.
The reader independently verifies144 structured fits and1149 intermediate
steps as well as the original probability/point/count checks. All four
joint-model controls select one population. At high depth,the state-plus-CNA
case selects two populations with120/120 globally matched assignments;
state-free whole MSE.00060328 improves pseudobulk.00363314 and matched
no-factor.00863210. The original minority criterion also passes. The complete
original mechanism gate is true. These are eight exposed simulation cases,
not replicated biological performance or a full-family success.

Low-depth true contrasts still select one population. The homogeneous joint
RNA variance cost remains. The fixed full290-case regression and unchanged
32/64 comparison families are therefore required;the eight-case pass is not
promotion. This is a model-estimation repair,not a change to false-split,
power,uncertainty or DNA-identification standards.

## Full regression execution

Writers23222366–23222373 cover the fixed290 cases in eight disjoint batches;
dependent reader23222374 requires every writer to succeed. All use nslab
account/QOS,CPU,1CPU8GB3h. The mathematical kernel,ordered minorizer and
fitter are byte-identical to the qualified mechanism version. The original
scientific comparison function is AST-identical to its frozen predecessor.
Independent reconstruction of every intermediate fit step is retained.
All290 predictions must lock before targets open. Full results pending.

[Frozen full design](../dev/ordered_amplitude_full_20261111/DESIGN.md),
[full execution identities](../dev/ordered_amplitude_full_20261111/SNAPSHOT.json),
[all job commands and reader dependency](../dev/ordered_amplitude_full_20261111/submissions.json),
[complete case population](../dev/ordered_amplitude_full_20261111/case_jobs.json).

## Full execution blocked by a traced optimizer bookkeeping defect

Writer23222370 fails at case147 under the unchanged nested-parent safeguard.
Its selected MAP558074.6077971692 is worse than an exact feasible parent
558054.8193714021,far beyond the.00882809 allowance. Replay23222406 reproduces
the failure and captures two ABNORMAL local optimizer attempts,including the
exact48-node parent evaluated at96 nodes. Discarding failed optimizer attempts
also discarded that better valid point. The other running writers are untouched.
Pending reader23222374 was cancelled because afterok can never succeed.
No full-study scientific results are available and no gate is weakened.

The bounded correction preserves evaluated parents as feasible incumbents after
at least one fixed-start solve succeeds. If a retained parent wins,it is marked
as lacking a new local-convergence certificate. All optimizer failures and the
original nested-parent safeguard remain. Mathematical count model,prior,starts,
scientific criteria and tolerances are unchanged. Preflight must independently
check the retained parent and pass the original eight-case gate before another
full run. This is an optimization bookkeeping repair,not proof of global
optimality or biological power.

[Exact failing fit and optimizer messages](../dev/ordered_parent_trace_20261111/run_23222406/results/FAILURE.json),
[replay source and original-source identity](../dev/ordered_parent_trace_20261111/trace.py),
[fixed scientific and computational contract](../dev/ordered_amplitude_retain_20261111/DESIGN.md).

## Incumbent-retention preflight passes;full corrected run resumes

Job23222482 completes0:0 in79s,scheduler peak73604KiB. Exact retained
parameters and independently reconstructed full PMF pass. Case147 now completes
with the original best parent atMAP558054.8193714021;its local-convergence
certificate is explicitly absent. The original eight-case mechanism gate
passes again. No optimizer tolerance,start or scientific criterion changed.

Eight corrected full writers23222513–23222520 and dependent reader23222521
run the unchanged290 cases and32/64 comparison families in new directories.
All use account/QOSnslab CPU,1CPU8GB3h. Summary23222569 depends on the reader
and uses1CPU8GB1h to preserve all failures,exact decision recounts,descriptive
means versus prior fitted candidates and every retained-incumbent status.
This summary adds no acceptance criteria. Full results pending.

[Exact preflight and retained-point status](../dev/ordered_amplitude_retain_20261111/preflight_results/REPORT.json),
[qualified current fitter](../dev/ordered_amplitude_retain_20261111/fitter.py),
[full corrected source identities](../dev/ordered_amplitude_retain_20261111/FULL_SNAPSHOT.json),
[all corrected full submissions](../dev/ordered_amplitude_retain_20261111/full_submissions.json),
[fixed descriptive summary](../dev/ordered_amplitude_retain_20261111/summarize.py).

## Complete corrected result:reject promotion

All corrected full writers,reader and summary complete0:0. Reader23222521 takes
13m57.02s,104448KiB process peak RSS;summary23222569 takes3.19s,44844KiB.
All5220 structured fits and41100 intermediate count/MAP updates independently
verify. Nine exact feasible-parent incumbents are retained with explicit
no-new-local-convergence status. Inputs,full PMFs,points,adaptive probes,
all64 exact decision recounts and matched baseline identities pass.

The scientific result fails26/32 original RNA and43/64 simultaneous comparisons,
below the preceding five-structure selector31/32+58/64. The eight-case mechanism
pass did not generalize even across these exposed simulations. Reject promotion;
keep the preceding selector and all frozen outcomes unchanged.

| Scenario/depth,48 seed units each | Ordered RNA whole MSE | Five-structure | Pseudobulk |
|---|---:|---:|---:|
| Clean/300 |.00431896|.00348635|.00413259|
| Clean/3000 |.000270759|.000202234|.00348458|
| State/300 |.00724314|.00652925|.0161965|
| State/3000 |.00135716|.000500701|.0155446|

At high depth,the joint route selectsK2 in41/48 seeds in both heterogeneous
conditions,versus48/48 in the original unrestricted model. Globally matched
assignments are5280/5760 clean and5503/5760 state;the latter falls below the
original5758/5760. Low-depth joint routes selectK1 in every48/48 seed in both
conditions. Homogeneous controls remainK1 in96/96 depth/seed cases,and both
state-only controls remainK1;this is not biological whole-specimen error control.

High-depth state-free whole MSE is.00147749 clean and.00157129 with state,
versus preceding.000142985/.000340140. The largest losses therefore cannot be
summarized as a harmless smoothing improvement. No genomic/uncertainty/power
claim changes. Biological n=0. Full production recovery2/5,false-split14/200
and selected correlated coverage54/1000 remain unchanged.

The next smallest decisive analysis uses cached parameters only:project each
existing unrestricted fit into the declared three-level family and compare
its same count-MAP target with the returned ordered fit. A better feasible
projected point proves missed optimization;absence of one does not prove
global optimality or identify the prior as the cause. Independently reconstruct
all resolved gains. This precedes any new fit,start,penalty or method change.

[Complete independent scientific reader](../dev/ordered_amplitude_retain_20261111/reader_23222521/results/REPORT.json),
[all failures,matched means and retained-parent statuses](../dev/ordered_amplitude_retain_20261111/run_23222569/results/SUMMARY.json).

## Cached sign-projection result:zero better points,not an optimization certificate

Job23223153 completes0:0 in85s,scheduler peak71736KiB. All290 input/lock
identities and5220 corresponding cached comparisons pass. The existing
median/sign projection supplies zero better feasible points beyond the
arithmetic allowance. Consequently,no independent gain probes are applicable.
This does not establish global optimization or prove that the model/prior
caused failure. The first projection itself uses the same initializer.

Source inspection gives a specific limitation:if the median lies between
feature values,the sign initializer assigns every feature to a positive or
negative level and none to neutral. Thus it initializes two levels rather
than fitting the declared three-level approximation. Whether this explains
any actual search failure remains unverified. The next bounded check is an
exact three-means projection on the profile's one-dimensional effect values,
with independent exact-rational qualification before cached target comparison.
There is no new state count,amplitude,penalty,threshold or fitted count model.

[Complete zero-witness result](../dev/ordered_projection_trace_20261111/run_23223153/results/REPORT.json),
[fixed sign-projection diagnostic contract](../dev/ordered_projection_trace_20261111/DESIGN.md).

## Exact projection proves search misses;minimal initializer correction in test

Qualification23223245 and cached trace23223246 complete0:0. Integer-lattice
cut selection agrees exactly with independent exhaustive Fraction assignments,
including ties,constant values,permutation and sign invariance. The fixed
counterexample's approximation SSE improves941→2 without a numerical cutoff.
The complete cache test finds1554/5220 better feasible points,all independently
verified using full count probabilities and the unchanged prior. This proves
search misses under the declared objective,not predictive superiority.

Among high-depth full joint96 fits,better feasible points occur in37/48 clean
and27/48 state cases. The sign projection can discard the neutral level;its
zero-witness result was not an optimization certificate. These are development optimization results;
they do not alter the failed26/32+43/64 performance result.

The smallest implemented correction replaces only the sign-based initial
projection with qualified exact three-means. The same PCA/two-means starts,
their number,parents,count kernel,prior,optimizer and tolerances remain. All
other fitter function ASTs and native-kernel/prior bytes match. No optimized
unrestricted initializer or additional start is added.

Preflight23223327 uses nslab CPU,1CPU8GB1h. It must pass the original eight-case
mechanism gate and improve both fixed first lexical high-depth full joint96
search witnesses outside that panel:case146(clean),case241(state). Independent
full PMFs assess the same objective. All ten bundles lock before original
mechanism truth opens. A failed preflight blocks expansion;no threshold or
criterion is weakened. Full predictive benefit remains unestablished.

[Exact projection qualification](../dev/ordered_three_means_20261111/qualification/REPORT.json),
[all1554 independently verified witnesses](../dev/ordered_three_means_20261111/run_23223246/results/REPORT.json),
[minimal corrected fitter](../dev/ordered_initialization_fix_20261111/fitter.py),
[fixed preflight and falsification](../dev/ordered_initialization_fix_20261111/DESIGN.md),
[exact source identities](../dev/ordered_initialization_fix_20261111/SNAPSHOT.json).

## Initializer-only correction passes preflight;full validation running

Job23223327 completes0:0. The original eight-case mechanism gate remains true.
Independent full-PMF target gains are14.942141 for case146 and20.489428 for
case241,against numerical allowances.0368404/.0710754. Thus the minimal
same-start projection change repairs both fixed optimization witnesses.
This is not yet a full predictive gain or a genomic membership result.

Writers23223372–23223379 and dependent reader23223380 now run the original
290 cases,unchanged32/64 comparison families,and full independent intermediate
step/probability/point checks into new output directories. All request nslab
account/QOS CPU1CPU8GB3h. No additional starts,count family,state count,prior,
threshold or optimizer tolerance is introduced. Preserve all previous
failures and the best preceding selector31/32+58/64. Do not promote from the
preflight. Full result pending.

[Completed preflight with both numerical comparisons](../dev/ordered_initialization_fix_20261111/preflight_results/REPORT.json),
[current full source identities](../dev/ordered_initialization_fix_20261111/FULL_SNAPSHOT.json),
[all current full commands and dependencies](../dev/ordered_initialization_fix_20261111/full_submissions.json).

## Live continuation and fixed full-summary contract

The next goal turn confirms all eight writers23223372–23223379 complete0:0;
reader23223380 is live on a compute node with no reported error. Summary23223425
requires its successful completion and uses nslab CPU1CPU8GB1h. The fixed
summary retains all32/64 decisions and failures,describes all eight strata and
compares pseudobulk,the prior five-structure selector,optional-state selection,
and the rejected ordered fit under identical case identities. It introduces
no new primary criterion or selected subset. Raw failed runs remain preserved.

Known-covariance uncertainty references remain conditional research results;
the independently identified reference/covariance and biological calibration
gaps are not repaired by this profile experiment. No new covariance/default
or uncertainty claim is added while this scientific gate is pending.

[Fixed descriptive summary source](../dev/ordered_initialization_fix_20261111/summarize.py),
[summary snapshot](../dev/ordered_initialization_fix_20261111/SUMMARY_SNAPSHOT.json),
[exact dependent summary submission](../dev/ordered_initialization_fix_20261111/summary_submission.json).

## Complete initializer-only result:high-depth repair,seven low-depth failures

All eight writers,reader23223380 andsummary23223425 complete0:0. The reader
takes12m07.12s,104024KiB process peak RSS;summary3.71s,60704KiB. All5220 fitted
model records and34802 intermediate steps independently verify. Four exact
feasible incumbents retain explicit missing-new-local-convergence status.
Every input,full PMF,point,adaptive check and64 exact decision recount passes.

The unchanged scientific family passes30/32 RNA and57/64 simultaneous tests,
up from26/32+43/64 before the initializer correction. All32 high-depth
comparisons pass. High-depth K2 selection recovers41/48→48/48 in both
heterogeneous conditions;globally matched assignments recover5280→5756/5760
clean and5503→5758/5760 with state. This is measured optimization/recovery
progress on exposed simulations,not genomic validation or a production change.

| High-depth whole-target mean MSE,48 seed units | Before initializer fix | After | Prior five-structure | Pseudobulk |
|---|---:|---:|---:|---:|
| Clean RNA |.000270759|.000274140|.000202234|.00348458|
| Clean state-free response |.00147749|.0000979238|.000142985|.00348458|
| State RNA |.00135716|.000498188|.000500701|.0155446|
| State-free response with state |.00157129|.000265158|.000340140|.00375644|

These means are descriptive paired-development summaries;there is no new
superiority test against the prior selector. The clean RNA mean does not
improve,despite the state-free gain. Do not combine favourable endpoints
into an all-task benefit.

Seven simultaneous failures remain,all atdepth300:clean whole RNA versus
pseudobulk24/48 wins;state minority RNA versus Beta19/48;clean minority
state-free versus global16/48 andBeta0/48;state minority state-free versus
pseudobulk0/48,global0/48 andBeta5/48. Clean low-depth K2 remains0/48;state
K2 becomes1/48 with global matched assignments4585/5760 versus4608/5760
before. That isolated extra distinction is not successful clone recovery.
Homogeneous96/96 and both state-only controls retainK1,but these do not
establish the full biological false-split target.

Reject general promotion under the unchanged complete gate. The prior
five-structure selector remains the best overall fitted regression31/32+58/64.
Public recovery2/5,continuous-state false splits14/200 and correlated selected
coverage54/1000 remain unchanged. Biological n=0. Do not add blind state/rank/
penalty variants to force a pass. A separate metadata-only feasibility gate
checks independent DNA/RNA development evidence before any genomic-response
learning;former evaluation specimens remain regression-only.

[Complete independent reader](../dev/ordered_initialization_fix_20261111/reader_23223380/results/REPORT.json),
[all failures,matched means and retained-parent records](../dev/ordered_initialization_fix_20261111/run_23223425/results/SUMMARY.json),
[independent DNA/RNA feasibility](DNA_RNA_RESPONSE_FEASIBILITY_2026-11-11.md).

The final checkpoint verifies all91 protected package/native/data/test identities,
all prior363 ledger values,14 immutable execution snapshots and114 metadata
input hashes. The public suite is not rerun because its source remains exact;
the historical521-Python/four-Rust results remain associated with that source.
No frozen result or scientific criterion changes. All ten full-study jobs are
terminalsuccess,with no remaining downstream dependency.
[Completed scheduler accounting](../dev/ordered_initialization_fix_20261111/accounting_complete.tsv),
[preservation checkpoint](../dev/ordered_initialization_fix_20261111/preservation_checkpoint_v2.json).
