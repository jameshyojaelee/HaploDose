# Count-native continuous-state safeguard: bounded performance attempt

The reference-only gene-split/weighted candidate recovered4/5 original cases but
failed the state-null error test (13/100). It remains rejected. This attempt
replaces its failed Gaussian log-ratio safeguard with an NB count observation
model that integrates a shared continuous per-cell regional log-rate effect.
Production clone inference is unchanged.

## Scientific contract and falsification

The estimand is a discrete RNA-profile population contrast conditional on a
Gaussian continuous state, not a DNA-defined clone. Both competing families use
the proposal's actual per-bin counts, expected rates and NB sizes. Null has one
location plus state SD; alternative has two locations, a mixture weight and
shared state SD. The existing BIC comparison charges two additional parameters.
This working-model criterion does not establish5% specimen error control.

The proposed correction is that count likelihood and shared latent state better
represent measurement uncertainty than the earlier Gaussian log-ratio model.
Observation-model misspecification, imperfect nuisance estimates and RNA-equivalent
biology remain competing explanations. Only independent corresponding DNA and
biological controls can establish the genomic interpretation. Biological n=0.

Gauss-Hermite quadrature uses24/48/96/192 points, with a declared convergence
and abstention rule. Null and alternative use multiple fixed starts including
the nested null. An independent numerical gate uses SciPy NB probabilities,
adaptive quadrature with explicit tail/rounding allowances, finite-difference
gradients and nested-null identity. Optimization is local; quadrature convergence
is not a formal integration enclosure or scientific confidence interval.

All original posterior.95, credible-core, gene/cell/effective support and30
recovery assertions remain. The safeguard acts only on otherwise accepted
proposals, so cannot rescue109's already rejected confirmation. First require
at least the retained4/5 recovery. Fewer passes reject the safeguard before a
large error/power expansion. Surviving recovery would still require fresh
complete-workflow error/power experiments matched to the mixed specimen regime.
Previous600-tumour/no-normal banks do not certify510-tumour/90-normal inputs.

## Execution status

[Frozen design](../dev/nb_state_recovery_20261109/DESIGN.md).
`nb_state_guard.py` implements the count likelihood; `nb_state_numerics.py` the
independent numerical gate; `nb_state_recovery.py` the complete original panel.
Both use immutable source_v3 and the new isolated noneditable candidate wheel
only after its full validation passes. Original unstarted jobs23137488/489 were
canceled after their build dependency failed a new list-input fixture. The fix
retains that test. Corrected numerical job23138005 precedes recovery23138006.
Both corrected jobs complete0:0. The four independent likelihood fixtures agree
to at most2.11e-15 with SciPy PMF/adaptive quadrature; finite-difference gradients,
quadrature-order checks and identical-component null identity pass. These are
computational checks, not a justification for the biological observation model.

## Complete recovery outcome and decision

| Seed | Pass every original assertion | Normal FPR | Tumour TPR | S1 groups | Assignment accuracy | Difference chromosomes |
|---|---|---:|---:|---:|---:|---|
|20261031|yes|.044444|.968627|2|.836032|chr12|
|20261107|no|.077778|.990196|1|.598020|none|
|20261108|yes|.088889|.992157|2|.756917|chr12|
|20261109|no|.022222|.980392|1|.602000|none|
|20261110|yes|.066667|.980392|2|.814000|chr12|

All cases retain one S2 tumour group; S3 also returns one group. The original
assignment denominator remains fitted tumour calls, assessed jointly with the
separate exactly-two-group requirement. None of the30 original assertions is
replaced. The safeguard retains the108/110 rescues but loses the previously
successful107. It recovers3/5 versus4/5 before guarding, failing the declared
prerequisite. **Reject this safeguard and do not launch the larger null/power
continuation.** Production remains2/5; the three production failures are not
resolved by this experiment.

All five tested safeguard comparisons resolve their BIC sign consistently at
orders24/48. Recorded advantages include9.6200, -7.9675,32.2506,-10.5697 and
8.8092, with changes at most about2.1e-9 across the paired quadrature orders.
Thus increasing integration precision cannot plausibly reverse the lost
comparison's negative8-unit advantage. Full records retain the source window,
cell identities, fitted parameters, objectives and convergence studies. The
109 case already rejects at weighted confirmation before this safeguard is
reached; a post-confirmation veto cannot fix it.

The scientific error/power tradeoff remains failed/unestablished: the new
method loses required recovery, and its false-split rate was **not estimated**.
No5% guarantee follows from the numerical tests or apparent simulation rescues.
Do not lower BIC penalties, credible-core support or the specimen error bound
to pass these cases. All failed attempts and frozen results remain intact.
[All30 criteria and verdicts](../dev/nb_state_recovery_20261109/recovery_v1/report.json),
[numerical references](../dev/nb_state_recovery_20261109/numerics_v1/report.json),
[execution records](../dev/nb_state_recovery_20261109).
