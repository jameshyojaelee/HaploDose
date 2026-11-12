# Discovery-predicted effects and gene membership do not rescue matched confirmation

Two implemented discovery-conditioned corrections fail the remaining109 core
gate. Refitting only the regional separation gives matched support.136048;
jointly fitting gene responsibilities and separation gives.024943. Both retain
an initial core and the other four paired examples, but neither gives109 a
matched core. Reject both as sufficient recovery fixes. No full recovery,
merge or larger false-split/power run follows. Production remains2/5.

## Scientific question and locked panel

Event discovery fits two regional RNA rate multipliers freely; confirmation
uses specimen beta*log(1.5). The [first design](../dev/predictive_gene_effect_20261109/DESIGN.md)
tests whether that effect-size mismatch loses confirmation evidence. The
[second design](../dev/predictive_gene_effect_20261109/MEMBERSHIP_DESIGN.md)
then tests whether aggregated-bin responsibilities disagree with the gene-level
observation model. Both designs precede their new method outcomes.

Estimand: regional relative RNA-rate separation between proposed populations,
conditional on the selected interval and nuisance fit. Working observation is
independent gene NB counts with the same reference-estimated sizes and expected
counts. Separation/membership must transport from discovery to confirmation.
Correlated expression state, gene composition and reference/noise mismatch are
alternatives. Independent corresponding DNA, suitable truth resolution and
complete-workflow null/power evidence establish biological validity, not a
conditional likelihood or an optimizer's success. Biological n=0 here.

Use ALL five already executed matched comparisons and their paired initial
comparisons in the fixed five-case count_gene experiment:20261031,20261107,
20261108,20261109,20261110. This is a mechanism panel, not five biological
specimens or a population recovery estimate. The complete180-entry original
trace and abstentions remain linked; other initial proposals are not silently
counted as new executions or reprocessed unnecessarily.

## Implementation and independent computation

[Research implementation](../benchmarks/hardening/predictive_gene_effect.py)
reconstructs discovery gene counts from the saved CSR proposal folds with
literal cell/gene IDs. NB sizes come from those same saved folds; expectations
are the saved proposal-derived expectations rescaled by the actual fold share.
Before proposal locking, the code accesses only proposal CSR entries and model
metadata, not full source CSR or held-out count entries. Input hashing still
records the full source files. Hashing is provenance, not use of test values.

First method retains the original q. Independently maximize each component's
responsibility-weighted NB likelihood on proposal observations only. A monotone
scalar score gives the unique finite rate when positive weighted counts support
it. Zero mass or numerically unresolved separation abstains; no rate floor.
Set the symmetric confirmation grid to[-delta,0,delta], where delta is the
absolute difference of fitted log rates. This is a conditional weighted-rate
estimate, not a fully new joint gene-discovery model or calibrated effect SE.

Second method jointly fits rates and responsibilities using the EXISTING
certified gene NB mixture optimizer, with unchanged initialization,mixing
bounds,score/observed-curvature and gain-error checks. Recompute component
proportions consistently with the new q and proposal-derived expected gene
counts in each bin. No new region search or starts are added. Local numerical
stationarity does not establish global optimality or biological identifiability.
Explicitly uncertified fits would be recorded as abstentions; all five actual
fits certify in this experiment. No convergence criterion is relaxed.

Each method locks ALL five parameter records before opening held-out gene
count arrays. Keep original genes/cells,expectations,priors,boundaries and the
same exact averaged beta-binomial likelihood-ratio model. Uniform-base matching
uses original full-count depth/eight strata, as in the failed.131252 baseline;
q remains integrated in the cell-shared likelihood. Minimum30 actual/effective
cells and every-bin.95/minimum5-bin cores remain. Probability mass cannot replace
actual cell support. The full workflow would still retain gene veto and merging.

Independent Fraction equal-exposure rates5/4 and7/4, Decimal80/100 scalar score
roots, component swap/permutation and confirmation-perturbation checks pass.
Maximum independent log-rate discrepancy2.22e-16. Bounds use operation count,
curvature,score residual and solver precision; they are numerical, not sampling
uncertainty. The second method additionally agrees with separately assembled
SciPy NB full-mixture log-PMF gain, discrepancy2.56e-13 within its derived
8.60e-9 allowance. Existing independent kernel/count-label-path checks remain.

All five original matching-weight calculations are exact before alternatives.
Replay controls evaluate original q/effect/uniform matching separately from new
q. Maximum posterior discrepancy across the ten controls is1.034e-25 in each
method;109's initial and matched controls are exactly equal. The tiny discrepancy
passes numerical consistency but is NOT literal exactness across the full panel.
Do not replace earlier failed literal or historical full-workflow controls.

## Actual109 result and all paired outcomes

| Model | Separation | Hard minority cells | Effective minority cells | Initial max support | Matched max support | Matched core |
|---|---:|---:|---:|---:|---:|---|
| Old q/step,uniform matching |.218899|31|31|.999899|.131252|None|
| Proposal weighted rates,old q |.232974|31|31|.999939|.136048|None|
| Joint gene q/rates |.253776|35|35|.999984|.024943|None|

The joint method increases matched total weight mass62 to70 and retains502
positive-weight cells, but stronger initial support does not transport into a
matched core. Membership polishing alone is therefore not a sufficient fix.
Both methods give a qualifying paired initial/matched core in the other four
locked examples. These are ten conditional comparisons per method, NOT4/5 full
clone recovery or an assessment of membership accuracy against DNA.

The continuation gate requires109's BOTH initial and matched core with the
unchanged support. Both fail. No historical full-control/candidate chain or
fresh null/power expansion is submitted. Existing production three recovery
failures, earlier4/5 state-null failure and simultaneous95% upper<=5% specimen
criterion remain unresolved. More signal in a conditional statistic cannot
substitute for that workflow or missing biological replication.

[Fixed-q results](../dev/predictive_gene_effect_20261109/panel_v1/report.json),
[joint results](../dev/predictive_gene_effect_20261109/joint_panel_v1/report.json),
[fixed-q locked discovery parameters](../dev/predictive_gene_effect_20261109/panel_v1/PROPOSAL_PARAMETERS_LOCKED.json),
[joint locked discovery parameters](../dev/predictive_gene_effect_20261109/joint_panel_v1/PROPOSAL_PARAMETERS_LOCKED.json).

## Execution and preservation

| Job | Role | Account/QOS; CPUs / RAM / requested time | Outcome |
|---|---|---|---|
|23151879|weighted-rate numerics|nslab;1 /4GB /1h|complete0:0,20s|
|23151880|five fixed-q proposal pairs|innovation;1 /8GB /1h|complete0:0,26s; gate fails|
|23152108|joint-mixture numerics|nslab;1 /4GB /1h|complete0:0,18s|
|23152109|five joint-q proposal pairs|innovation;1 /8GB /1h|complete0:0,24s; gate fails|

CPU compute node006,neutral cwd,existing read-only noneditable validated wheel
and native identity. Peak process RSS352668KiB. Commands,allocations,backend,
versions,resources,dependencies and snapshots are saved in
[the task directory](../dev/predictive_gene_effect_20261109/). Source_v1 runs
the first pair; source_v2 is unsubmitted; source_v3 adds the explicitly planned
uncertified-optimizer abstention before executing the joint pair. There is no
failed optimizer run being silently retried. Shared installations,production,
native code,original tests,criteria and frozen results remain unchanged.
No new full-package/GPU test or biological calibration is claimed.
