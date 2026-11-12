# Reference-estimated gene noise: one additional recovery, no error guarantee

The reference-estimated, weighted gene-count candidate passes **3/5** original
recovery cases, versus current production **2/5**, without the known simulation
noise parameter. It newly recovers20261110. The two other failed cases remain.
Fresh paired pilot simulations retain10/10 two-profile recovery but show1/10
expression-state false splits in both methods. The candidate is not promoted.
Production therefore remains2/5; frozen results and all criteria are unchanged.

## Implemented estimator and scientific assumptions

[reference_gene_noise.py](../benchmarks/hardening/reference_gene_noise.py) estimates
one common gene-level NB dispersion from declared fit-reference RNA only.
It uses no query cell, held-out reference, confirmation outcome or true-noise
constant. Reference groups have separate gene means. For each group, with counts
Y_cg and exposure L_c, estimate squared rate using cross-cell products:

    r_g²_hat = ((sum_c Y_cg)² - sum_c Y_cg²) /
               ((sum_c L_c)² - sum_c L_c²).

Under independent NB observations with mean L_c*r_g, the cross-cell term avoids
the same-cell variance contribution. Pooling groups and genes gives:

    phi_hat = sum_tg [sum_c Y_cg(Y_cg-1) - r_tg²_hat * sum_c L_c²] /
              sum_tg [r_tg²_hat * sum_c L_c²].

The numerator/denominator moments have the intended expectations under known
exposures and common dispersion. Their finite-sample ratio is not unbiased.
Actual exposures are used-gene library totals, so compositional dependence and
reference cell-state heterogeneity remain working-model limitations. Positive
finite fitted dispersion is required; there is no oracle fallback or outcome-
selected clipping. All training row IDs and group moments are saved.

The estimand is common residual RNA count dispersion under that working model,
followed by relative regional RNA effects and recoverable populations. It is
not DNA dosage or clone identity. Gene-specific noise, latent reference state,
inferred exposures and reference mismatch can violate transport to query cells.
Independent technical sampling measurements would distinguish count noise from
biological heterogeneity; corresponding DNA remains necessary for genomic claims.

[gene_thinning_hook.py](../benchmarks/hardening/gene_thinning_hook.py) now accepts
a callable estimated size alongside its retained historical scalar size path.
The selected size is saved with actual gene splits; gene-to-bin count conservation
and veto-half correspondence remain required. The scalar historical oracle arm
was not rerun this turn; its frozen results/source remain untouched. No production
package, native source, dependency or shared installation changed.

## Full original recovery experiment

The [frozen design](../dev/reference_gene_noise_20261109/DESIGN.md) defines two
arms: estimated gene splits with hard pooling, and those splits with the existing
both-stage membership weighting plus effective-unit matched gene veto. The
unchanged original six assertions score all five cases, not a selected subset.

| Source case | Production | Estimated gene/hard | Estimated gene/weighted |
|---|---|---|---|
|20261031|pass|fail|pass|
|20261107|pass|fail|pass|
|20261108|fail|fail|fail|
|20261109|fail|fail|fail|
|20261110|fail|fail|pass|
|Total|2/5|0/5|3/5|

Both new arms estimate identical phi values per case, approximately.1642–.1756
(sizes5.70–6.09). They do not reuse the known-size10 oracle. In newly recovered
20261110, the weighted arm has normal FPR6/90, tumour detection500/510, exactly
two S1 tumour groups, one S2 group, assignment accuracy.812 and chr12-only
distinctions. All original thresholds pass. These are exposed development cases,
not five independent biological specimens or a general recovery guarantee.

The gain is an interaction of estimated gene splitting, weighting and the matched
gene veto; this experiment does not isolate a causal benefit of dispersion
estimation alone. The0/5 hard arm is retained prominently. All parameters,
criteria, gene halves and actual weights are under
[hard_v1](../dev/reference_gene_noise_20261109/hard_v1) and
[both_v1](../dev/reference_gene_noise_20261109/both_v1).

## Fresh full-workflow false-split and power pilot

The [pilot design](../dev/reference_gene_noise_20261109/PILOT_DESIGN.md) freezes
the weighted candidate before new simulations. It uses the existing primary
gene-count generator with fresh phase2, replicates0–9, in three conditions:
one CNA population under the fitted law; one with correlated expression state;
and two CNA populations with the same nuisance. Production and candidate run
complete workflows on identical raw inputs, including reference fitting,
initialization, all proposals, gene evidence, merging, splitting and stopping.
The generator's noise is never passed to the candidate.

| Condition | Production | Candidate |
|---|---:|---:|
|Matched-null false splits|0/10|0/10|
|Expression-state-null false splits|1/10|1/10|
|Exact chr12 two-profile recovery|10/10|10/10|
|Two-clone mean all-cell membership accuracy|.72250|.72733|

Ten simulated workflow replicates per condition are the replication units, not
6000 cells. Biological n=0. Membership is scored against generated labels;
the pilot adds no real same-cell DNA validation. Mean matched-null accuracy is
.99867 in both; state-null accuracy is.99000 production versus.96667 candidate,
so the negative membership result also remains visible.

The frozen pilot rule allocates alpha.05/(2 nulls ×2 looks)=.0125 per condition/
look. Exact binomial excess-error tails would reject at3 or more false splits
out of10. One error gives p=.401263, so this pilot does **not** statistically
reject a5% population rate. It also does not establish error control: the
98.75% upper bounds are.354805 for0/10 and.490636 for1/10.

The originally contemplated total100 panel could not certify the target for
this candidate even with no further state-null errors: with1/100, the exact
upper bound at this confidence level exceeds5%. Equivalently,
P[Binomial(100,.05)<=1] is greater than.0125. Stop that planned certification
continuation for futility. This is not proof that the population rate exceeds5%,
and the horizon/criterion is not enlarged or loosened after seeing outcomes.
The candidate remains a partial recovery improvement, not a validated release.

[All pilot performance](../dev/reference_gene_noise_20261109/pilot_summary_v1/performance.json),
[summary](../dev/reference_gene_noise_20261109/pilot_summary_v1/report.json),
[pilot implementation](../benchmarks/hardening/reference_gene_noise_pilot.py) and
[exact performance collector](../benchmarks/hardening/reference_gene_noise_summary.py)
retain all60 paired-arm workflows, raw inputs, reference fits, labels and outputs.
No error or power case was dropped. There is no stronger-baseline profile win or
spatial claim from this clone experiment.

## Execution and remaining work

Jobs23131383/384 run the original hard/weighted recovery arms. Jobs23131434–439
run six pilot shards;23131491 collects results. All nine complete on CPU nodes,
using innovation/nslab, immutable source and the validated noneditable wheel.
The eight inference jobs request4CPU/16GB/2h; the collector1CPU/4GB/2h. No arrays
or interactive QOS. Source_v1 executes original recovery; source_v2 pilot;
source_v3 collector. The candidate estimator and gene-split implementation stay
unchanged across those snapshots.

[Exact job commands](../dev/reference_gene_noise_20261109/jobs.json),
[scheduler accounting](../dev/reference_gene_noise_20261109/accounting.tsv) and
[process resources](../dev/reference_gene_noise_20261109/resource_summary.json)
record timings, memory and exits. The previous403-Python/4-Rust production check
is historical; the full suite was not rerun for research-only changes. Actual
original recovery and pilot workflows were executed instead.

Continue performance development on the remaining20261108/109 information/model
failures and the state-null distinction. Retain this non-oracle3/5 candidate for
comparisons. It cannot replace production or be described as completing the
three-regression repair, meeting5% control or establishing SOTA.
