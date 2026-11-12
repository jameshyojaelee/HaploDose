# Five original recoveries, followed by a failed false-split pilot

The experimental depth-conditioned/gene-score workflow passes **all5 original
recovery cases**, including the3 production failures. The unchanged production
package still passes2/5. The experimental candidate then falsely splits4/10
continuous-state null simulations, versus production0/10, and fails the
prespecified error rule. All three procedures recover10/10 two-population
simulations. Keep the recovery milestone and error failure together: this
candidate is rejected for promotion. No5% specimen guarantee or biological
clone validation follows. No frozen results or scientific criteria changed.

## What was changed, and why

The [preceding full workflow](DEPTH_CONDITIONED_WORKFLOW_2026-11-09.md)
retained109's eight-bin confirmation core but rejected its matched-margin gene
veto,46/73 signs,p=.017208 at unchanged.01. Switching to a passing hard pool
would not be justified. Instead implement a gene statistic for the **same**
frozen proposal responsibilities, depth strata, common-support cells and raw
gene observation model used in the new regional confirmation.

[Prespecified design](../dev/conditioned_gene_score_20261109/DESIGN.md),
[full-workflow design](../dev/conditioned_gene_score_20261109/FULL_DESIGN.md),
[gene statistic](../benchmarks/hardening/conditioned_gene_score.py) and
[experimental integration](../benchmarks/hardening/gene_shared_confirmation.py).
No production/native/original-test/dependency or shared-environment edits.

Estimand: gene RNA contrast along a selected mixture after gene-by-depth-stratum
abundance nuisance. For each gene/stratum, fit its null NB log rate on the
confirmation observations, inside the gene test. With p=mu/(mu+r), efficient
contrast score is sum[(q-qbar)*(y-(y+r)*p)], qbar=sum(r*p*q)/sum(r*p).
Information is sum[r*p*(q-qbar)^2]. Proposal q, cells, genes and strata remain
locked before confirmation; the gene-test nuisance fit never feeds discovery.

Count a score sign only beyond its NB root-error/operation-rounding bound.
Unresolved signs remain in the denominator, not concordant. No within-stratum
q variation provides no contrast information; zero-count null-rate boundaries
give zero score/information. Preserve informative count>=10, min8 genes and
nominal one-sided sign p<.01. Insufficient genes now abstain rather than
automatically pass. Preserve every-bin.95/min5-bin/min30-cell support, callability,
reference-only quintile noise, weighted confirmation/merge and stopping.
Only the active supported event proposal uses the new gene statistic; merge
gene checks retain their existing implementation. Both old/new decisions are
recorded. No hard-pool substitution or merge bypass.

Identifying assumptions: working independent NB genes, adequate depth strata,
reference expectations and transport of reference-only gene noise. Regional
expression state, gene dependence and reference mismatch can produce the same
gene concordance. Nominal sign p is **not** a specimen error probability. The
adaptive hypothesis family is every searched proposal and veto in the workflow.
The fresh pilot below tests this failure directly; independently corresponding
DNA at appropriate resolution is still required for a genomic interpretation.

## Independent computation and original recovery

Exact rational reference score36/77 and information41/77 agree with computation
to5.55e-17 and0, respectively, within derived4.00e-11 and2.56e-13 budgets.
Independent SciPy full NB log-PMF finite differences at4 step sizes pass their
third-derivative truncation and rounding bounds. Exact binomial tails, component
swap, gene/cell permutations, pure-depth/no-q-information controls, zero-count
boundaries and insufficient-gene abstention pass. These are numerical checks,
not sampling calibration. [Numerical report](../dev/conditioned_gene_score_20261109/numerics_v1/report.json).

Before full execution, lock all5 qualifying previously executed matched cores
across all5 original cases, retaining all180 source trace entries and previous
decisions. All5 new gene checks pass, including all4 previously accepted cores.
For109 the frozen eight-bin core has56/78 supporting signs,p=7.47408e-5.
[Complete bounded panel](../dev/conditioned_gene_score_20261109/panel_v1/report.json).

The original-runtime correction and exact historical tests are preserved. The
new-source control reproduces every original array/table/15 fold and all30
scientific criteria exactly before the full candidate runs. Captured runtime003,
the read-only validated noneditable wheel and neutral compute directories are
used; no tolerance is weakened. Earlier failed controls remain unchanged.
[Strict new-source control](../dev/conditioned_gene_score_20261109/full_control_v2/report.json).

| Original seed | Pass every original assertion | Normal FPR /90 | Tumour TPR /510 | S1 groups | Assignment accuracy, original denominator | Difference |
|---|---|---:|---:|---:|---:|---|
|20261031|yes|.044444|.968627|2|.836032 /494|chr12|
|20261107|yes|.066667|.988235|2|.738095 /504|chr12|
|20261108|yes|.088889|.992157|2|.756917 /506|chr12|
|20261109|yes|.022222|.970588|2|.800000 /495|chr12|
|20261110|yes|.066667|.980392|2|.814000 /500|chr12|

All retain one S2 and S3 group. Original assignment criterion> .65 is interpreted
jointly with exactly2 groups and chr12-only differences; its fitted-call
denominator is not relabeled all-cell accuracy.109's regional minimum support
remains.999649655, mean.999766770; only the population-aligned gene test changes
its rejection. Weighted merging remains active and the contrast survives it.
[All original criteria](../dev/conditioned_gene_score_20261109/full_candidate_v2/report.json).

## Fresh whole-workflow error and simultaneous power

Freeze all candidate/transitive method hashes before fresh phase6 outcomes.
[Pilot design](../dev/conditioned_gene_score_20261109/PILOT_DESIGN.md),
[driver](../benchmarks/hardening/conditioned_gene_pilot.py),
[prespecified summary](../benchmarks/hardening/conditioned_gene_pilot_summary.py).
Three methods share each raw dataset: production, depth_conditioned with old
gene veto, depth_gene_score. Exact reference-noise and raw gene/proposal-fold
identity is checked between research arms. Truth fields are removed from model
metadata; generation truth is opened by the scorer after each fit.

10 independent base seeds generate30 paired condition datasets and90 complete
workflows, including nuisance fitting, initialization, proposal selection,
confirmation, gene checks, refinement, merge and stopping. Conditions: one CNA
population under the model; same population with continuous regional RNA state;
two CNA populations with the same state. Both nulls and methods are paired,
not independent biological specimens. Power effect/prevalence remain fixed
delta=.6*log(1.5),40%; no new biological sample is added.

| Method | Matched null errors /10 | State-null errors /10 | Exact chr12 recovery /10 | Mean membership accuracy, all600 cells |
|---|---:|---:|---:|---:|
|Production|0|0|10|.737833|
|Depth-conditioned, old gene veto|0|3|10|.739833|
|Depth-conditioned, new gene score|0|4|10|.739833|

The primary candidate's error family remains2 nulls x2 looks (pilot10/final100),
one-sided exact alpha.0125 per null/look, Bonferroni simultaneous95%. Its state
excess-tail p=.00102850 is below1/80; simultaneous upper bound is.771840.
Even with no subsequent errors its final100 upper would be.108748>.05.
Reject and stop the planned certification continuation; do not enlarge the
horizon after outcomes. All power fits finish regardless.10/10 power has an
ordinary descriptive exact95% interval[.691503,1], not guaranteed100% recovery.
No additional candidate-only power gains occur in these paired10 seeds.
No method has zero reported groups. Zero matched-null errors give an upper
.354805 at this pilot size, not5% certification. Production's pilot0/10 state
errors also do not establish biological error control.
[All error bounds, membership, power and abstention](../dev/conditioned_gene_score_20261109/pilot_summary_v3/performance.json),
[decision](../dev/conditioned_gene_score_20261109/pilot_summary_v3/report.json).

Three state failures are shared by both research methods; one is introduced by
the new gene score. The nuisance/model change already supplies most unsupported
distinctions, and the new statistic does not resolve state specificity. This
is direct counterevidence to promoting five-case recovery alone. Do not tune
sign alpha, credible support, merge thresholds or the specimen criterion.

## Execution and preservation

| Job | Role | Account/QOS; CPU /GB /time | Outcome |
|---|---|---|---|
|23153776|independent numerical checks|nslab;1 /4 /1h|complete0:0,24s|
|23153777|locked all5-core panel|innovation;1 /8 /1h|complete0:0,5s|
|23153846|unchanged strict full control|innovation;4 /16 /3h|complete0:0,90s|
|23153847|full5/5 recovery candidate|innovation;4 /32 /3h|complete0:0,143s|
|23154218|matched null0..4,3 arms|innovation;4 /16 /2h|complete0:0,135s|
|23154219|matched null5..9,3 arms|nslab;4 /16 /2h|complete0:0,172s|
|23154220|state null0..4,3 arms|innovation;4 /16 /2h|complete0:0,253s|
|23154221|state null5..9,3 arms|nslab;4 /16 /2h|complete0:0,196s|
|23154222|power0..4,3 arms|innovation;4 /16 /2h|complete0:0,219s|
|23154223|power5..9,3 arms|nslab;4 /16 /2h|complete0:0,207s|
|23154224|dependent exact summary|nslab;1 /4 /1h|complete0:0,8s|

All CPU compute nodes, no arrays/interactive QOS. Immutable source_v1 numerical/
bounded work,source_v2 full controls/candidate,source_v3 adds pilot drivers only;
all previous method files are exact betweenv2/v3. Exact submissions, snapshots,
versions, actual backend, `/usr/bin/time` resources and failures remain in the
[task directory](../dev/conditioned_gene_score_20261109/).
Full package/GPU suites are not rerun for research-only changes. Biological
false-split calibration and membership validation are not run. The experiment
fixes the locked recovery regression in a research mode, while failing the
simultaneous scientific performance/error objective. Production remains unchanged.
