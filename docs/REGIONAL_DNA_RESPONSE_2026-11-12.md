# Pool-before-fit prediction of measured relative DNA

This is a concrete profile-performance experiment, not another correctness audit
or a claim that the full HaploDose goal is achieved. It tests whether pooling
nearby RNA measurements before fitting their response improves prediction of
processed relative DNA. Original recovery, whole-specimen false-split/power,
uncertainty, independent-evaluation and fixed-competitor requirements remain.

The subsequent [ordered mechanism trace and matched ablation](REGIONAL_RESPONSE_CAUSE_2026-11-12.md)
preserve that outcome. Removing RNA absent from its selected training reference
fails the separately fixed mean-and-median diagnostic gate; stop that route.
No new method or rescue tuning is promoted from the accounting stages.

## Completed conditional-null result and final decision

All **99/99** fixed permutations and summary complete 0:0. Writer jobs
23234026/23234027/23234028 take 9m03s/9m01s/8m58s; timed process peaks are
545,164/544,328/544,704 KiB, with scheduler samples 495,556/462,852/470,872 KiB.
Summary23234029 takes16 Slurm seconds, timed15.49s, peak84,480KiB. All four
120-identity before/after manifests match. All nine observed arrays reproduce
bitwise in each writer; each null executes the exact same fitting block and
passes its leakage, independent OLS and primary scalar-loss checks. That is
9,792 OLS probes including the three observed replays and 99 trials, plus
106,029 checked primary unit losses in the trials. Every full prediction lock
and array hash verifies in the summary.

There are 97 movable and 22 singleton units, and
46,167,057,931,450,682,815,272,714,240,000,000 possible within-stratum pairings.
The observed minimum mean advantage is0.00242382438, with numerical envelope
1.82e-13. **0/99** null statistics reach its conservative lower envelope.
Randomization **p=0.01**; the exact-binomial-method one-sided95% Monte Carlo
tail upper is approximately **0.02980667377**. Independent80-digit Decimal
zero-exceedance calculation differs from the numerical SciPy quantile by three
binary64 ULPs, far from the fixed0.05 continuation threshold.

This supports conditional association under the declared exchangeability null
for the fixed, exposed workflow. It does not replay historical model selection,
test the composite claim of superiority to each baseline, establish biological
false-split control, or replace independent validation. The 99 permutations
do not add biological specimens. The eight-comparator paired-effect assessment
below still fails. **Final decision: preserve the point gain and conditional
association, withhold promotion and broad superiority, and do not tune a rescue
on these outcomes.** The consumed50 units remain excluded from development.

[Final null report](../dev/regional_dna_null_20261112/summary_23234029/results/REPORT.json),
[exact fitting/scoring workflow](../dev/regional_dna_null_20261112/workflow.py),
[first complete writer](../dev/regional_dna_null_20261112/run_23234026/results/COMPLETE.json),
[summary resources](../dev/regional_dna_null_20261112/summary_23234029/resources.txt).

## Paired-effect result: superiority is withheld

The follow-up **23234048 completes 0:0** in 7 Slurm seconds, timed 6.98 seconds,
peak RSS 200,496 KiB. Its 117 before/after source identities match. The unequal-
cluster rational fixture and all 320,000 resampled effect calculations pass;
maximum scalar discrepancy is 1.11e-16. There are 119 units in 19 lineage
clusters, not 320,000 biological replicates.

The candidate improves the previous gene-specific predictor in **39/119** units,
ties in 18, and has median paired gain **−0.0001285511**. Fold 2 regresses by
mean **−0.0008270730**. These adverse outcomes remain in the complete result.
Both simultaneous working lower bounds against that predictor cross zero:
**−0.0021122223** under unit resampling and **−0.0049634542** under lineage
resampling. Only the two naive comparators pass both working bounds; 8/16
individual bounds are positive. Strict superiority over the fixed family fails.

The original mean-error continuation result is preserved, but it cannot justify
promotion, a broadly improved typical-unit claim, or an independent performance
claim. The fixed 99-permutation experiment continues to its declared horizon;
even a passing conditional association result cannot override these bounds.
No alternate pooling, favorable lineage, confidence method or sample size is
selected to rescue this candidate. No consumed test values are reopened.

These working bounds condition on fitted models and do not cover historical
selection, shared training/reference estimation or unverified biological
dependence. They are evidence limiting the claim, not newly calibrated package
uncertainty. Original selected-interval coverage remains 54/1,000.

[Complete paired-effect result](../dev/regional_dna_uncertainty_20261112/run_23234048/results/REPORT.json),
[saved indices and resampled effects](../dev/regional_dna_uncertainty_20261112/run_23234048/results/REPRODUCTION.npz),
[exact implementation](../dev/regional_dna_uncertainty_20261112/run.py),
[resources](../dev/regional_dna_uncertainty_20261112/run_23234048/resources.txt).

## Completed development result: a small profile gain

Job **23233997 completes 0:0** in 32 Slurm seconds. The scientific process takes
24.29 seconds, with peak RSS 483,820 KiB; scheduler sampled peak is 406,608 KiB.
All **109** before/after source identities match, including 91 protected source,
native, data and test identities. All 11 original and five new fixtures pass.
All five preceding prediction arrays replay bitwise. Held-out DNA perturbations,
96 independent SciPy OLS references, all integer read sums, and all 2,142
primary/secondary unit-loss calculations pass.

All 119 development units, 17,751 gene loci and 2,361 occupied 1 Mb regions
remain. Every method uses the same final truth support and operators.

| Fixed method | Mean 1 Mb MSE | Mean arm MSE |
|---|---:|---:|
| Previous gene-level reference | 0.295591624 | 0.069980558 |
| Previous gene-level naive | 0.355782184 | 0.092532319 |
| Previous gene-level global | 0.294620317 | 0.069425683 |
| Previous gene-specific response | 0.282871134 | 0.067590349 |
| Previous nested blend | 0.283507505 | 0.067666890 |
| Regional reference | 0.295591624 | 0.069980558 |
| Regional naive | 0.367223029 | 0.076140963 |
| Regional global | 0.294210181 | 0.068699715 |
| **Regional-specific response** | **0.280447309** | **0.064896870** |

The candidate passes the frozen eight-comparator point-mean continuation gate.
Its smallest primary mean advantage is 0.00242382438 against the previous
gene-specific response, versus an arithmetic budget of 1.58e-13. That budget
describes numerical calculation, not sampling uncertainty. Pooling alone is
not uniformly beneficial: regional naive prediction is worse than its previous
gene-level counterpart. The result does not establish genomic causation,
same-cell performance or independent superiority.

[Full result](../dev/regional_dna_response_20261112/run_23233997/results/REPORT.json),
[all unit losses](../dev/regional_dna_response_20261112/run_23233997/results/DONOR_LOSSES.tsv),
[independent checks](../dev/regional_dna_response_20261112/run_23233997/results/INDEPENDENT_CHECKS.json),
[original fixtures](../dev/regional_dna_response_20261112/run_23233997/ORIGINAL_FIXTURES.json),
[new fixtures](../dev/regional_dna_response_20261112/run_23233997/results/FIXTURES.json),
[resources](../dev/regional_dna_response_20261112/run_23233997/resources.txt).

## Frozen follow-ups: association and paired-effect uncertainty

The [complete-workflow conditional-null design](../dev/regional_dna_null_20261112/DESIGN.md)
retains all nine methods and repeats every fitting step under 99 fixed whole-RNA
row permutations within lineage/library/DNA-assay strata. No test unit or method
is dropped. The exact observed fitting block is compiled from its frozen AST;
every writer first reproduces all nine observed prediction arrays bitwise.
Each draw retains independent OLS, leakage and primary unit-loss checks, full
predictions, parameters and a mapping journal. Historical method development
is not replayed, and conditional exchangeability remains an assumption.

The single statistic is the minimum candidate mean advantage across all eight
baselines. Roundoff overlap counts conservatively as a null exceedance.
The fixed 99-draw horizon yields a randomization p-value and exact one-sided
95% Monte Carlo tail upper bound. That is association evidence for this fixed
workflow, not a biological false-split bound or eight separate superiority tests.

Three nslab writers **23234026–23234028** each request CPU 1 CPU, 4 GB, 2 hours;
afterok summary **23234029** requests 1 CPU, 2 GB, 30 minutes. All use v6&amd,
120 immutable identities and the original isolated environment/runtime. At the
first checkpoint, all three nine-array observed replays pass and 51/99 draws
are complete; no aggregate null outcome is inspected. No arrays are submitted.
[Executed workflow](../dev/regional_dna_null_20261112/workflow.py),
[source inventory](../dev/regional_dna_null_20261112/SNAPSHOT.json).

A [separate paired-effect protocol](../dev/regional_dna_uncertainty_20261112/DESIGN.md)
was written after observing means but before examining per-unit wins or bounds.
It reports medians, wins and fold effects, then 20,000 unit and 20,000 lineage
resamples for 16 Bonferroni-protected one-sided basic working bounds. Unequal
clusters retain equal-unit weighting. All 320,000 resampled effects are checked
independently with scalar sums. No confidence type, horizon or subgroup is
selected from its output. Shared fitted-training/reference uncertainty and
historical selection are not covered; biological coverage remains unverified.
Any nonpositive or unresolved bound withholds a strict superiority claim,
regardless of the association-null outcome. **23234048** requests nslab CPU
1 CPU, 4 GB, 30 minutes, with 117 immutable identities. Its outcome is pending
at this checkpoint. These follow-ups do not reopen the 50 consumed test units.

## Motivation and implemented intervention

The earlier [gene-response experiment](DNA_RNA_RESPONSE_DEVELOPMENT_2026-11-11.md)
beats its simple development baselines but its
[consumed50-unit test](DNA_RNA_RESERVED_2026-11-11.md) fails strict superiority.
Those50 units remain regression-only and are not used for development here.
The current study uses the original119 development units and three fixed folds.

Gene-specific RNA noise can obscure a shared regional DNA signal. Sum the
retained RNA counts within the already fixed1Mb regions, then fit the unchanged
training-only reference and response machinery at that level. No new width,
penalty,reference group,noise constant or tuning search is introduced. RNA is
normalized by the same retained-read library total. Regional training targets
are mean processed linear relative DNA over the fixed gene loci. Predict each
held-out unit with its RNA and lineage/library metadata only. Broadcast regional
predictions to genes before applying the unchanged common-support scorer.

For independent Poisson genes with fixed proportions and one regional multiplier,
the sum is sufficient. This motivates the intervention; it does not assert a
Poisson law for biological RNA or count-model uncertainty for this regression.
Highly expressed state genes, dosage compensation, internal DNA breakpoints and
lineage/processing confounding can defeat this pooling assumption.

The candidate is regional-gene OLS. Eight fixed comparators remain: the preceding
five gene-level methods(reference,naive,global,gene,nested blend), and three
matched regional baselines(reference,naive,global). The latter receive exactly
the same regional support,reference construction and projection. A gain over
the former alone cannot establish that regional OLS adds information.
This is not a substitute for the original single-cell pseudobulk benchmark.

## Scientific contract and falsification

The primary estimand is equal-unit mean within-unit1Mb MSE for measured linear
relative DNA at the same finite gene loci. Arm MSE is secondary, obtained from
the same fitted1Mb predictions rather than a separately optimized arm model.
PatientID-labelled units, not genes or windows, define the repeated-unit level;
their independent biological provenance remains unverified.

The working model is a training-only linear conditional predictor. Identifying
assumptions are no donor/reference leakage and transport of the RNA/DNA relation
across the fixed development folds. Regional transcription,culture and technical
processing can predict DNA without establishing genomic causation. Independent
corresponding DNA with suitable assay resolution and biological replication,
and same-cell truth for membership, remain the distinguishing measurements.

The frozen continuation gate requires lower primary mean error than **all eight**
comparators beyond numerical uncertainty. Ties,unavailable outputs and unresolved
arithmetic cannot pass. Failure rejects this candidate without a follow-on
width,subgroup or parameter search. A pass permits only separately frozen
workflow-null and independent validation; no public promotion follows here.

## Exact implementation, checks and execution

[Frozen design](../dev/regional_dna_response_20261112/DESIGN.md),
[executable study and independent calculations](../dev/regional_dna_response_20261112/run.py),
[launcher](../dev/regional_dna_response_20261112/run.sbatch),
[immutable source inventory](../dev/regional_dna_response_20261112/SNAPSHOT.json).

The study retains the original11 numerical/identity/support/leakage fixtures,
adds hand-counted pooling,constant-region,missing-DNA and permutation fixtures,
and requires exact conservation of integer RNA counts. All five previous
prediction arrays must replay bitwise. Held-out DNA replacement withNaN must
leave fits and predictions identical. Predictions lock before held-out scoring.
Ninety-six independently formed full-design SciPy OLS checks and scalar
reconstruction of every primary/secondary loss use operation/conditioning error
budgets rather than a universal epsilon. Numerical correctness cannot establish
biological identification or superiority.

Job **23233997**: nslab account/QOS,CPU partition,1CPU4GB2h,originalv6&amd
constraint,existing isolated read-only Conda/glibc runtime. No array,install,
public edit,shared-environment write,frozen overwrite or running-job change.
Initial static AST and shell syntax checks pass;109 immutable identities include
91 protected package/native/data/test identities. At this initial checkpoint,
the job is running on`ne1dc6-007.nygenome.org`; no final result is claimed.

Submission:

```sh
sbatch --parsable --output=/gpfs/commons/home/jameslee/Haplodose/dev/regional_dna_response_20261112/slurm-%j.out dev/regional_dna_response_20261112/run.sbatch /gpfs/commons/home/jameslee/Haplodose/dev/regional_dna_response_20261112
```

The four full-DNA producers23233127–23233130 continue independently. Their
native reader23233877 remains dependent on successful completion; none of these
running sources or allocations is changed. No RNA/DNA truth claim follows from
their live status. The original production recovery2/5,state false splits14/200,
power195/200 and selected correlated coverage54/1000 remain unchanged.
