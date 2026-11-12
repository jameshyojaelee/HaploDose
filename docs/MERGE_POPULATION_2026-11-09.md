# Merge populations: separate membership dilution from depth matching

The diagnostic asks where the accepted chr12 distinction in recovery case20261108
is lost: when proposal-half memberships define the merge populations, or when
depth balancing changes their composition. It changes no inference method.
Original recovery remains2/5; no candidate correction or false-split/power result
has been established by this study.

## Contract and decisive comparison

The estimand is the conditional test-half regional RNA-rate contrast between the
actual assigned cell pools. The observation model is the executed depth-balanced
pooled count model. Pure genomic clone effects are a different estimand. Selected
membership mixing, depth-associated clone composition and RNA-model mismatch can
each reduce the observed contrast; correlated expression state remains a biological
alternative. Same-cell DNA at adequate sensitivity would distinguish genomic
identity, whereas this diagnostic only localizes the computation/population step.

The [design](../dev/merge_population_20261109/DESIGN.md) retains all five original
cases and all actual merge pairs. It captures exact test counts, rows, labels,
proposal/full responsibilities, incoming doses, actual weights, exposures,
first/second depth moments, squared-weight variance counts and returned segments.
Truth is joined only after inference for population composition summaries.

Compare unweighted and depth-matched first10Mb chr12 centered RNA log-rate means
on identical bins. These means use an explicit descriptive half-count pseudocount
and genome median centering; they are not inferential probabilities or DNA-copy
effects. Report clone1/clone2 fractions and total pool weights, retaining other
cells through the denominators. Empty support is explicitly untestable.

A loss only after matching implicates depth-associated composition; weak
separation already beforehand implicates membership. Neither change would instead
direct attention to the observation/core model. These patterns are mechanism
evidence, not automatic justification for a new default or biological claim.

## Independent gates

The observer must preserve every public output and all30 original recovery
criteria exactly. It must capture nonempty comparisons. The independent audit
reconstructs pooled counts and squared-weight variance counts with scalar
compensated sums, plus exposure and depth moments. Numerical guards scale with
reduction length, binary64 epsilon and summed term magnitude. For positive count
sums the bound is32*epsilon*n*max(sum,1); signed depth sums use sum of absolute
terms. The factor32 covers product/reduction ordering conservatively; these are
checked floating-error allowances, not statistical tolerances or confidence
bounds. Empty audits cannot pass the corrected gate.

No posterior, core length, group size, seed, comparator, frozen result or
false-split criterion changes. A correction, if justified later, still needs
complete-workflow false-split and detection-power validation.

## Attempts and execution

| Job | Resources | Observed outcome |
|---|---|---|
|23130578|innovation CPU,4CPU/16GB/2h|Completed0:0 in1m31s; exact public/30-criterion parity, but captured0 comparisons. Mechanism investigation invalid. Scheduler batch MaxRSS993472K.|
|23130579|nslab CPU,1CPU/4GB/2h|Completed0:0 in5s, but audited0 pairs/0 scalar checks. Vacuous software pass is not accepted evidence.|
|23130609|innovation CPU,4CPU/16GB/2h|Completed0:0 in1m13s;4 captured comparisons; exact public/30-criterion parity.|
|23130610|nslab CPU,1CPU/4GB/2h|Completed0:0 in6s;4 pairs and6696 scalar count/variance checks, exposures and depth moments pass.|
|23130718|nslab CPU,1CPU/4GB/2h|Completed0:0 in5s; saved full-MAP versus proposal-MAP population follow-up, all4 pairs retained.|

The original observer expected the original matched-difference function directly
under the merge frame. The existing membership-stage runner wraps that function
even in its hard control. The frame condition silently excluded actual calls.
The [amendment](../dev/merge_population_20261109/OBSERVER_AMENDMENT.md)
walks actual ancestors to the original merge method, still requiring the immediate
caller to be the original matched-difference function. Both trace and audit now
require nonempty observations. The original zero-capture outputs remain intact.

Code: [merge_population.py](../benchmarks/hardening/merge_population.py).
The executed zero-capture source is source_v2; source_v1 is an unused intermediate.
The replacement is source_v4; source_v3 is an unused intermediate before adding
explicit amendment provenance and transitive input/output verification.
[Exact commands and outputs](../dev/merge_population_20261109/jobs.json)
and [task artifacts](../dev/merge_population_20261109) preserve both attempts.
All jobs use the previously validated noneditable wheel/isolated environment,
immutable source, no arrays and no interactive QOS. No production code changes.

## Observed population changes and interpretation

The corrected trace contains four merge comparisons: one each in the two passing
cases20261031/20261107, and two in failed case20261108, all specimenS1. The other
two failed cases do not reach this merge-pair branch; they remain initial
confirmation failures. Captured results preserve the observed supported segments
in the two passing cases and empty segments for both20261108 merge comparisons.

The original [audit_v2](../dev/merge_population_20261109/audit_v2/report.json)
and the subsequent [audit_v3](../dev/merge_population_20261109/audit_v3/report.json)
pass. Their full population/contrast records are in
[comparisons.json](../dev/merge_population_20261109/audit_v3/comparisons.json).
The [follow-up design](../dev/merge_population_20261109/MEMBERSHIP_FOLLOWUP.md)
was written after observing the depth-matching result and before opening the
full-MAP comparison. It preserves all four pairs, no new inference and no new
selection. Source_v5 executes only the additional saved-array audit.

| Case/comparison | Clone1 composition separation: full-MAP / proposal-MAP / depth-matched | Absolute centered chr12 RNA contrast: same order |
|---|---|---|
|20261031|.658661 / .443180 / .445752|.300158 / .123848 / .115816|
|20261107|.677126 / .540849 / .554572|.356333 / .173979 / .163303|
|20261108 first|.466668 / .412516 / .419899|.298835 / .115512 / .099969|
|20261108 second|.576598 / .460773 / .482080|.279437 / .087552 / .075427|

Composition separation is the absolute difference of clone1 fractions between
the pools, with all cells retained in denominators; it is not a classification
accuracy or inferential statistic. Contrast values are RNA log-rate differences.
Full-count assignments use test evidence. Their larger contrast can reflect both
membership information and selection on those same observations. It cannot be
reported as independently confirmed signal or as a causal effect of only better
membership. These are four dependent simulation comparisons, biological n=0.

[Observed pattern, not a guarantee] The main recorded contrast reduction precedes
depth balancing. Matching slightly improves true-clone composition separation in
all four comparisons, while modestly reducing the descriptive RNA contrast.
Thus loss of clone-composition separation *during depth matching* does not explain
these failed merges. The data do not rule out other consequences of matching,
including effective information, nuisance effects or the final core model.

Retaining full-count memberships is not a justified automatic repair. That
joint-data mechanism previously reached5/5 recovery and failed the fresh
[false-split/power experiment](JOINT_FRESH_TRADEOFF_2026-11-08.md). This new trace
does not erase that negative evidence or supply a new calibrated procedure.
Next, a justified correction must model uncertain memberships/observation scale
through the actual merge statistic and pass the complete workflow under matched
and misspecified nulls alongside power. A conditional fixed-pool model probe can
reject a mechanism before such an expensive workflow test; favorable full-data
contrasts cannot substitute for it.

## Provenance and limits

All successful replacement input/output manifest file hashes were independently
rechecked, including both actual20261108 comparisons. Exact public output parity
and all30 original criteria include the three failures; recovery stays2/5.
The full production suite was not rerun because this is observer/audit-only code.
No fresh null or power experiment was run and no new biological specimen was
added. Five jobs exited0, but the first two are scientifically rejected zero-
capture attempts, not five successful evidence checks.

See [scheduler accounting](../dev/merge_population_20261109/accounting.tsv) and
[process resources](../dev/merge_population_20261109/resource_summary.json).
Process peak RSS is1307752/187404/1311404/196452/197376KiB for the five jobs in
table order. Process and scheduler memory reports differ and are retained
separately. Frozen inputs/results, production methods, native source and tests
remain preserved; no new dependency, installation or default was introduced.
