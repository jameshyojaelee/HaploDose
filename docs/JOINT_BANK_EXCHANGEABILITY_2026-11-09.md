# Frozen false-split excess remains consistent with finite-bank imbalance

An exact conditional diagnostic does not reject exchangeability of the frozen
matched-law calibration/evaluation scores after Holm adjustment across16 checks.
This is consistent with finite calibration-bank imbalance. It does not establish
that mechanism, prove the observation model correct, or turn the failed5%
false-split/power criterion into a pass. No correction or operating point is
promoted and no additional inference draws are generated.

## Question, observation and assumptions

The estimand is a conditional upper-tail probability for the observed number of
accepted evaluation specimens, given the399 frozen scores:199 calibration and
200 evaluation per arm/null. Under exchangeable scores, every199-element
calibration subset has equal probability. Both joint arms, both matched NB/state
null pairs and all four existing gates are retained, giving16 post hoc diagnostic
hypotheses with exact-rational Holm adjustment. The scientific family153 and all
frozen scientific confidence bounds remain unchanged.

This is a cause diagnostic after observing the earlier failures, not independent
confirmation. Algorithm/role mismatch or different generating laws could violate
exchangeability; finite sampling alone can also produce apparent tail imbalance.
The earlier source-traced workflow audit addresses implementation correspondence;
this conditional distributional check cannot prove those assumptions statistically.
Neither interpretation identifies genomic clones or substitutes for biological
replication. Biological n=0; cells/seeds/score allocations do not add specimens.

## Exact calculation and actual findings

For each observed score, the original rank rule uses all calibration scores
greater than or equal to it. The calculation counts every calibration allocation
using descending tied-score blocks, retaining zero-score exclusion. It uses exact
integer combination counts; no floating tolerance determines a result.64 small
fixtures agree with independent exhaustive subset enumeration. Unique positive
score examples also satisfy an independent hypergeometric counting identity.

All16 original event counts reproduce. Selected entries below are the two
primary state diagnostics; the complete16 results are retained in exact_v1.

| Arm, original gate.05 | Original state false splits | Nominal diagnostic tail | Holm16 tail |
|---|---:|---:|---:|
| Joint proposal-membership control |19/200|.03921697|.58825461|
| Joint full-membership candidate |20/200|.02694659|.43114540|

Both strict.005 state rules retain3/200 events, nominal tail.12499763,Holm1.
Matched-NB rules retain0/200 events and diagnostic tail1. No diagnostic rejects
exchangeability under its declared family. Non-rejection is not evidence of
equivalence, correct coverage, adequate power or confirmed finite-bank causation.
In particular, these p-values are not estimates or upper bounds for the specimen
false-split rate. The original simultaneous rate-bound failure stands.

## Independent actual-score reference

Audit23033146 verifies all32 input score-file identities and reconstructs all16
events using sorted calibration scores and bisect_left. For every positive
observed count s, the actual pooled boundary at m=s+h is strictly positive and
not crossed by ties, where h=200*alpha-1. Only under that checked condition,
K>=s is equivalent to at most h calibration labels among the first m ranks.
Independent integer hypergeometric sums reproduce every actual diagnostic tail
exactly as fractions. Zero observed counts have tail1. All16 allocation masses
and rational Holm results agree exactly, independently of the tied-block dynamic
program. A failed boundary condition would invalidate this reference, not justify
ignoring ties. No such failure occurs.

## Decision and remaining evidence

Do not change the three locked recovery assertions, loosen the scientific gate,
add draws to the frozen experiment or declare a discovered calibration bug.
Current evidence does not justify another score/threshold correction from this
panel. A future calibration procedure needs a separate prespecified assessment
including calibration uncertainty and simultaneous recovery; it cannot rewrite
these results or establish biological control from simulations alone.

The RNA/state identification failure remains: independently justified nuisance
constraints or suitable same-cell DNA measurements are necessary for stronger
genomic interpretation. Available raw libraries in the existing scONE cell-line
candidate may support a development measurement-chain investigation after fixing
its role and assay contract; processed brain matrices do not establish its raw
DNA chain. No new cohort is admitted here. DEFND/OneK1K correspondence/resource
questions and the biological null-panel gap remain unresolved.

## Execution and changes

Exact diagnostic23032687 completes13s,processpeak89876KiB; audit23033146
completes8s,processpeak89692KiB. Both exit0 on nslab CPU nodes, using the current
isolated validated wheel read-only. No package, test, scientific criterion,
environment or frozen result changes. Full suite not rerun for diagnostic-only
scripts; latest evidence remains365 Python+4 Rust pass,1 GPU-only skip.
`git diff --check` passes. Prior-wheelGPU23024875 still queued, no GPU pass claimed.

Added scripts:`benchmarks/hardening/joint_bank_exchangeability.py` and
`joint_bank_exchangeability_audit.py`. Protocols, immutable source_v2 and
audit_source_v1, exact probability distributions, input identities, full reports,
commands and actual environment/resource records are under
`dev/joint_bank_exchangeability_20261109/{exact_v1,audit_v1}`. The unused source_v1
is preserved. All original failure results and decisions remain unchanged.
