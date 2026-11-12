# Soft dispersion moments do not repair recovery or state-driven splitting

Hard membership changes the population target of the shared dispersion moment,
even under a correctly specified NB mixture. A responsibility-weighted version
preserves that ideal-model identity, but it does not rescue the three recovery
failures or improve the fresh false-split/power experiment. Keep the counterexample
and research implementation; do not promote this method as a recovery or
calibration correction. No production source, original test, frozen output or
scientific criterion changes.

## The assumption tested

The estimand is conditional RNA count dispersion in
Var(Y|k)=mu_k+phi*mu_k². With known correct component means and posterior
responsibilities, iterated expectation gives

    E[sum_k R_k(Y)*((Y-mu_k)^2-mu_k)] /
    E[sum_k R_k(Y)*mu_k^2] = phi.

Replacing R by a MAP label conditions residuals on the observed count's selected
component. It need not preserve this target. This is a statistical moment
distinction, not proof that the package's heuristic violates a claimed exact EM
likelihood maximization. Fitted means, estimated memberships, normalization and
model misspecification remain additional assumptions in an actual RNA workflow.

The strongest competing explanation for apparent genomic subdivisions remains
correlated transcriptional state. A per-bin dispersion correction cannot by
itself distinguish that state from a CNA-shaped mean shift. Corresponding DNA
and independently justified state/reference information would be distinguishing
biological measurements. Simulations below establish no biological validity.

## Independent known-mixture counterexample

[Initial design](../dev/dispersion_membership_20261109/DESIGN.md) uses NB shape5,
phi0.20, mean pairs(4,6), (4,12), (4,30) and prior probabilities0.1,0.5,0.9 for
the lower-mean component: nine population calculations, biological n=0.
Exposure is independently fixed; counts are not normalized by their own total.

The initial calculation subtracts the observed count in its Poisson correction.
Production subtracts the fitted mean. Both have the same soft-weight population
identity, but different selected-residual moments, so
[the matched follow-up](../dev/dispersion_membership_20261109/MATCHED_DESIGN.md)
repeats all nine cases with the exact production subtraction. The initial run
and its successful audit remain unchanged; it is not substituted for the
production-matched result.

For the production-matched calculation:

| Means | Lower-component prior | True/soft phi | Hard-MAP moment phi |
|---|---:|---:|---:|
|4,6|0.1|0.20|0.194444|
|4,6|0.5|0.20|0.079775|
|4,6|0.9|0.20|0.252167|
|4,12|0.1|0.20|0.192549|
|4,12|0.5|0.20|0.123331|
|4,12|0.9|0.20|0.189379|
|4,30|0.1|0.20|0.198518|
|4,30|0.5|0.20|0.184949|
|4,30|0.9|0.20|0.178466|

The direction is not uniformly downward. These are population targets, not a
finite-sample bias estimate. The calculation deliberately excludes clipping and
median shrinkage to isolate the membership step.

SciPy probability sums use NB factorial-moment identities to bound omitted
tails, doubling truncation until the requested ratio-tail accuracy is below
1e-12, then doubling again for refinement. This is a numerical request, not a
scientific threshold. Separate exact Fraction recurrences verify all18 matched
ratios and4617 MAP decisions without SciPy PMFs. Soft moments agree with0.20
within the certified tail/numerical allowance. The initial alternative
subtraction has its own18-ratio/4617-decision exact audit.

Code: [moment calculation](../benchmarks/hardening/dispersion_membership.py),
[rational audit](../benchmarks/hardening/dispersion_membership_audit.py).

## Matched implementation and original recovery

[Research hook](../benchmarks/hardening/dispersion_weighted.py) changes only
the dispersion moment's latent weights. For the same depth/programme-adjusted
means it accumulates sum R*((Y-m)^2-m) and sum R*m². It retains the original
divide fallback0.05, clipping[0.001,10] and half-per-bin/half-median shrinkage.
Kept-component weights are not renormalized: doing so would change this selected
latent-component target. The ideal identity is not a guarantee when components
or means have been changed during fitting.

Every update is checked through an independent polynomial expansion using
E_R[m], E_R[m²] and row mass, with an absolute-term cancellation bound. A local
MAP update at the same R/doses is recorded separately before restoring phi and
applying the soft update. It confirms that the hook executes; it is not confused
with the paired complete-run control.

[Frozen recovery design](../dev/dispersion_membership_20261109/RECOVERY_DESIGN.md)
and [runner](../benchmarks/hardening/dispersion_recovery.py) retain all five seeds,
raw-count inputs, truth-free fitting metadata, original assertions and six
scientific criteria. Current MAP and soft arms both pass2/5. Seeds20261108,
20261109 and20261110 still fail. The numerical moment distinction does not explain
enough of the lost membership information to repair these cases.

[Independent recovery audit](../benchmarks/hardening/dispersion_recovery_audit.py)
rejoins original identities, recounts all60 criteria, verifies thin mode, confirms
that original test outcomes match those criteria, and reconstructs each saved
soft dispersion update. All checks pass; the scientific recovery gate fails.

## Fresh complete-binned-workflow false-split and power test

[Prespecified design](../dev/dispersion_membership_20261109/TRADEOFF_DESIGN.md)
uses200 fresh paired simulations in each of three frozen same_nuisance_v2
conditions,600 cells/240 bins and20-bin chromosome blocks. The two arms differ
only in the moment update. Initialization, merging, splitting, stopping,
terminal validation and all thresholds remain unchanged. No joint-membership
research method is substituted for the current thin procedure.

All1200 fits complete. New seed keys are[20261109,501,condition_index,replicate];
all observed/nuisance/truth/reset hashes match between arms. Truth is not supplied
to fitting. The reference and initial dispersion are known simulation inputs.
This exercises the complete bin search but omits public raw-input preprocessing,
gene and allele channels. It is not a biological or public-workflow guarantee.

| Endpoint | MAP | Soft moment | Simultaneous bound, MAP / soft |
|---|---:|---:|---:|
| Matched-NB null, extra group |0/200|0/200|Upper0.03964 /0.03964|
| Correlated-state null, extra group |190/200|191/200|Upper0.98713 /0.98948|
| Two populations, multigroup detection |200/200|200/200|Lower0.96036 /0.96036|
| Two populations, exactly two groups |0/200|0/200|Lower0 /0|
| Exactly two and all-cell assignment>0.65 |0/200|0/200|Lower0 /0|

An extra group is counted even when final pair evidence is absent; missing
support cannot erase a grouping error. Every cell remains in assignment
denominators. Reference-like/extra groups receive no extra one-to-one credit.
The original0.65 assignment cutoff is an explicitly retained diagnostic, not
a universal biological power floor. Multigroup detection alone is plainly
insufficient: neither arm ever returns exactly two groups in these alternatives.

The ten new one-sided exact bounds are simultaneous with all153 historical
endpoints at Bonferroni alpha0.05/163. Their previous153-family bounds are retained
in a new record; old outputs are not overwritten. The null upper-bound criterion
remains at most0.05. Paired changed-decision counts are descriptive; no new
uncorrected superiority screen is introduced. Neither arm passes the state-null
criterion, and no general calibration or recovery claim is made.
Historical exposed analyses remain regression/development evidence; expanded
family accounting does not convert them into confirmatory validation.

[Collector](../benchmarks/hardening/dispersion_tradeoff_summary.py) requires every
completed chunk and paired input identity. [Independent audit](../benchmarks/hardening/dispersion_tradeoff_audit.py)
exhaustively verifies all1200 one-to-one assignment matchings, reconstructs
endpoint decisions and paired changes, and checks all163 bounds by60/80-digit
binomial-polynomial inversion. It verifies all153 historical event counts and
prior bound values unchanged. Computational checks pass; biological n remains0.

## Execution and next decision

Twenty-one ordinary CPU jobs across innovation/nslab complete with exit0; no
arrays or interactive QOS are used. Isolated wheel/source/backend identities,
resource records, dependencies and all immutable output paths are in
[jobs.json](../dev/dispersion_membership_20261109/jobs.json). The original73
package/native/test/README files remain unchanged. Production suite was not
rerun for this research-only method; the isolated inference studies and
independent audits above were executed on compute nodes.

Do not promote the weighted moment as a recovery or error-control fix. Retain
its ideal-model identity as a useful distinction for future variance modelling,
without claiming fitted conditional uncertainty is calibrated. Membership
formation and misspecified correlated state remain consequential. The three
original recovery failures, biological validation, supported clone membership
and competitive accuracy goals remain open. More simulations of this same
failed correction are not the next scientific step.
