# Research gene veto depends on arbitrary weight units

The aligned-pool research veto fails weight-unit invariance. A fixed exact
counterexample and three captured real-pool comparisons change veto decisions
when every pooling weight is multiplied by the same positive factor. The
regional conditional log-rate and variance remain identical. This blocks
promotion of the current research candidate, regardless of its5/5 exposed
recovery result or the pending fresh experiment. Production is unchanged.

## Estimand, law and decisive check

The target is the direction of relative gene abundance in a fixed weighted
population. Scaling W by f preserves that population and its relative abundance.
For fixed weights under the conditional Poisson sampling approximation,
S=sum(wY), Q=sum(w²Y), L=sum(w exposure), c=sum(w²)/sum(w).
Changing weight units gives(fS,f²Q,fL,fc). Installed pooled_lograte uses
log((S+.5c)/L) and (Q+.5c²)/(S+.5c)², which preserve their values.
This algebra does not establish NB coverage, independent weights or selection
validity. It identifies a required unit correspondence within that working model.

The research gene veto instead applies a fixed .5 pseudocount and eligibility
S_a+S_b>=10 to fractional weighted counts. Its numerical results therefore
depend on the absolute weight units. The strongest alternative is real loss of
counting information; it is excluded here because no count or cell changes and
the conditional mean/variance summaries are invariant under the same rescaling.

The locked24-gene fixture has8 regional genes, four with counts(20,10),four with
(10,20), and16 outside genes(20,20). At unit weights,4/8 regional signs agree,
exact binomial tail163/256, so the inherited .01 veto rejects. Multiplying all
weights by1/8 leaves the profile and conditional variance unchanged, but no
gene meets the unscaled10-count eligibility. The inherited insufficient-evidence
branch passes. That is a concrete unit-dependent support decision.

## All captured matched pools, fixed scales

The design fixes dyadic factors1/8,1/2,1,2,8 before execution. It inspects every
captured matched pool in the frozen original-five aligned candidate. Three
decision changes occur; other captured decisions remain recorded.

| Original case / call | Scale1/8 n/k, decision | Scale1 n/k, decision |
|---|---|---|
|20261107 / merge23|8/7, reject|54/40, pass|
|20261108 / confirmation30|10/9, reject|58/46, pass|
|20261109 / merge61|8/7, reject|60/45, pass|

The table reports n eligible genes and k concordant signs. These repeated
observations are development mechanisms, not independent specimens or an error
rate estimate. No preferred scale is selected. The original5/5 result and
running fresh candidate retain their original weights, code and thresholds.

For every pool, squared-weight count sums and c reconstruct the regional
conditional summaries. All five scales give exactly zero observed log-rate and
variance discrepancy in these probes. Independent scalar formulas agree within
a64-operation binary64 allowance; exact dyadic scaling avoids introducing a
new approximate scientific threshold. Rational centered-odds comparisons and
integer combination tails independently reproduce every veto decision. Ties
and zero-centered odds retain the existing rule, rather than being rounded away.

## Execution and decision

One CPU job23120542, nslab1CPU/4GB/1h, completes0:0 in19s. Design, source,
input hashes, results, exact tails and execution/resource records reside under
`dev/gene_veto_units_20261109`. Implementation:
`benchmarks/hardening/gene_veto_units.py`. The check reads only existing evidence;
it changes no public code/test, frozen output, running experiment or scientific
criterion. Full package tests are not rerun for the diagnostic. Biological n=0.

Do not promote the current weighted veto. A future correction must derive
eligibility and pseudocount units from the weighted sampling model, recover
the raw-count rule for binary unit weights, preserve common-weight rescaling,
and independently re-evaluate recovery and complete-workflow error/power.
Choosing a scale because it rescues a case is prohibited. Correct units alone
will still not calibrate selected gene signs, address fitted-membership/reference
uncertainty, or distinguish transcriptional states from CNA-defined clones.
