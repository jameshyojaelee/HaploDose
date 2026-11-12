# Regional-response gain: mechanism accounting and rejected correction

The lower development mean error is real, but the proposed removal of RNA
unsupported by its fitting reference does not pass the fixed mechanism test.
No correction is promoted. The earlier failed paired-effect superiority
assessment remains binding. This investigation explains a consequential
measurement distinction without relabelling it as improved genomic accuracy.

## Scientific question and fixed evidence

The estimand is the equal-unit mean change in 1 Mb processed-relative-DNA MSE,
conditional on the existing three-fold fits, references and RNA range guards.
The observation model is the existing working linear RNA/DNA response. There
are 119 development PatientID-labelled units and 2,361 regions. Those labels
do not certify independent specimens, same-cell truth or calibrated DNA states.
The 50 consumed test units remain excluded from development.

Pooling can change gene weighting, fitted slopes, support decisions and
nonnegative projection together. Sampling zeros, lineage, transcriptional state
and genomic dosage can each produce RNA absent from a selected reference.
Corresponding independent DNA and replicated RNA would distinguish these
explanations; arithmetic agreement cannot.

## Attempt 1: ordered accounting on frozen fits

[Design](../dev/regional_response_cause_20261112/DESIGN.md) and
[implementation](../dev/regional_response_cause_20261112/run.py) were frozen
before execution. Job **23234126 completes 0:0**, nslab/nslab, CPU1/4 GB/1 hour,
on ne1dc6-010. Slurm elapsed 10 seconds; timed process 8.68 seconds, maximum
RSS 286,772 KiB. The scheduler's 920 KiB sample misses the short process peak.
All **117** before/after identities match, including inherited protected source.

Both endpoint prediction arrays and fitted parameter digests reproduce exactly.
The hand-derived rational fixture, pooled-contrast identities, scalar stage
losses and per-unit telescoping checks pass. Largest pooled identity discrepancy
is 3.64e-12 against a maximum input-size/magnitude operation budget of 1.12e-8.
The endpoint accounting budget is 2.66e-12. These are arithmetic bounds, not
sampling intervals or scientific effect thresholds.

| Fixed ordered stage | Mean 1 Mb MSE |
|---|---:|
| Original guarded nonnegative gene predictions | 0.282871134 |
| Remove gene nonnegative projection | 0.283527490 |
| Remove gene RNA range guards | 9,377.369680 |
| Regional slopes with equal gene weights | 39,720.350596 |
| Training-reference RNA weights | 142.438964 |
| Include RNA absent from selected training reference | 282.049763 |
| Restore actual regional RNA guard | 0.281934834 |
| Restore actual regional nonnegative projection | 0.280447309 |

The endpoint gain remains 0.00242382438. Large opposing changes outside training
support make this ordered decomposition highly sensitive to interactions. It
does not provide a unique causal attribution, show that the gain is merely noise
averaging, or license any intermediate stage as a method. No subgroup or stage
is selected using its DNA performance.

One algebraic distinction survives that limitation: reference-weighted gene
contrasts do not reconstruct a pooled RNA contrast when query RNA comes from
genes with zero abundance in the selected training reference. This contribution
is separate from ordinary reweighting and can pass a regional guard even when
the gene learner abstains.

[All outcomes and numerical checks](../dev/regional_response_cause_20261112/run_23234126/results/REPORT.json),
[per-unit accounting](../dev/regional_response_cause_20261112/run_23234126/results/UNIT_ACCOUNTING.tsv),
[locked stage profiles](../dev/regional_response_cause_20261112/run_23234126/results/PROFILES_LOCKED.json),
[resources](../dev/regional_response_cause_20261112/run_23234126/resources.txt).

## Attempt 2: isolate that contribution while retaining actual guards

[Separately frozen design](../dev/regional_unseen_ablation_20261112/DESIGN.md)
holds the actual regional fits, references, full retained-read library
denominators, guard masks, nonnegative projection and all 119 units fixed.
Only the unsupported RNA contribution is removed. It does not refit a
reduced-feature estimator or renormalise after discarding reads. Counterfactual
profiles lock before held-out DNA scoring. The RNA-selected example's DNA
annotation was also placed after that lock before submission; no executed
failure was hidden or criterion relaxed.

The prespecified diagnostic gate requires numerically resolved positive **mean
and median** improvement. A failure stops this mechanism route; even a pass
would only justify a separately frozen intervention, not overturn superiority.

Job **23234253 completes 0:0**, nslab/nslab CPU1/4 GB/1 hour, ne1dc6-007.
Slurm elapsed 8 seconds; timed process 6.61 seconds, peak RSS 140,788 KiB.
The batch RSS sample is zero and does not estimate the process peak. All
**127** before/after identities and independent scalar losses pass.

| Fixed diagnostic | Result |
|---|---:|
| Actual regional MSE | 0.280447309416 |
| Fixed-guard counterfactual MSE | 0.280443221901 |
| Mean improvement | 0.000004087515 |
| Median improvement | **0** |
| Resolved unit wins / losses / ties or unresolved | 54 / 43 / 22 |
| Fold 1 mean improvement | 0.000006335435 |
| Fold 2 mean improvement | 0.000008658824 |
| Fold 3 mean improvement | **−0.000003612195** |

The mean/median arithmetic budgets are 2.66e-12 and 8.41e-12 respectively.
**The diagnostic gate fails. Reject this proposed recovery mechanism and stop
the route.** Do not choose a convenient subgroup, change a guard or retain the
tiny favourable mean while ignoring the zero median and adverse fold.

[Complete result](../dev/regional_unseen_ablation_20261112/run_23234253/results/REPORT.json),
[all unit losses](../dev/regional_unseen_ablation_20261112/run_23234253/results/UNIT_LOSSES.tsv),
[locked profiles](../dev/regional_unseen_ablation_20261112/run_23234253/results/PROFILES_LOCKED.json),
[resources](../dev/regional_unseen_ablation_20261112/run_23234253/resources.txt).

## Primary count trace selected without query DNA

The largest unrepresented RNA fraction in a region admitted by its RNA guard
selects ACH-002693 / PT-DyxJsT, RNAProfileID PR-anqrOi, fold 1. The gene is
PRAC1, ENSG00000159182, HGNC:30591, NCBI 84366; GRCh38 source interval
chr17:[48,721,718,48,722,518), region chr17:48, arm 17q. No biological gene
function or CNA is inferred from its identity.

The exact source RNA column is 46348, labelled
`PRAC1 (ENSG00000159182)`. Its RNA-SeQC read abundance is 15,680 of retained
118,181,987 reads. These are not verified independent molecules. The selected
Peripheral Nervous System / stranded reference contains ACH-001302, ACH-001338,
ACH-001354 and ACH-002066, with RNAProfileIDs PR-CAxmUa, PR-BoNJN3, PR-hOYnpE
and PR-0dHyNy. Every one has source gene count zero.

Only after the RNA profile lock, the already selected example is annotated:
processed gene/region DNA 0.834206235; training regional DNA 1.560675848;
actual prediction 1.416110327; fixed-guard counterfactual 1.372235411. The DNA
column is 19059, `PRAC1 (84366)`. These uncertain processed measurements neither
identify absolute copy state nor establish biological validation. An illustrative
local improvement cannot rescue the failed complete-panel mechanism gate.

## Decisions and remaining limits

The [original regional-response result](REGIONAL_DNA_RESPONSE_2026-11-12.md)
retains its mean gain, conditional association and failed superiority bounds.
No package/native/default/test source or frozen result changes. Original
production recovery remains 2/5, state-null false splits 14/200, selected
correlated coverage 54/1,000. The fixed competitor family and all criteria remain.
These studies resolve neither those failures nor missing biological replication.

The four original full-library DNA mappings remain untouched. Bounded stderr
inspection finds reference chrM extraction warnings in trim9 (read
SRR7214625.1655972_1655972/1, chrM 1) and trim0
(SRR7214625.2956130_2956130/1, chrM 2). chrM is outside the primary autosomal
profile; other-contig counts remain recorded. These warnings do not establish
an autosomal defect, their cause or clean completion. Preserve them pending
terminal native reconstruction; do not alter the running analyses.
