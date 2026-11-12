# Reference-only regional nuisance: strong conditional gains, failed safeguard

The normal-reference predictor improves all 64 prespecified state-challenge
comparisons, with 32/32 strict seed wins in every comparison. It fails six of
12 signal-preservation safeguards. The complete gate therefore fails: no real-data
expansion, production promotion, changed cutoff or claim of genomic validation.
Both nslab CPU jobs complete successfully. This investigates finding 5 of the
[scientific review](scientific_code_review_2026-11-10.md); it does not close it.

The scientific goal remains accurate genomic gain/loss profiles across cellular,
spatial and bulk RNA, with reliable CNA-defined populations and memberships when
independent measurements support those interpretations. Narrower current claims
describe the evidence available; they do not replace that goal.

## What changed and why

This is a distinct research estimator, confined to
[its development directory](../dev/within_region_nuisance_20261110/DESIGN.md).
No package, native code, original test or dependency changes accompany it.
The preceding failed scONE response and clone-prevalence candidates remain closed.

Current programme construction removes arm means from nuisance loadings. That
also leaves broad transcriptional state in RNA dosage. Here, independent normal
RNA learns a map from gene contrasts *within* regions to the state's regional
means. In the ideal log-intensity model, within-region centering removes both a
constant regional dosage effect and the global depth/normalizer. The predicted
state retains regional means, unlike a high-pass target. The matched full-expression
predictor is an ablation that can absorb dosage.

The estimand is the planted regional natural-log RNA response, centered over
regions, with an identity dosage-response bridge specified by simulation. Gene
counts are Poisson, optionally with heterogeneous independent gene Gamma rate
variation. Count integrality is known by construction; these are not deposited or
reconstructed measurements. The main identifying assumption is that informative
within-region state markers and their relationship to regional means transfer
from normals to queries. A coherent state without markers and a reversed
normal/query relationship are explicit competing explanations. Gene-specific
dosage response or breakpoints within a region would also defeat the ideal
projection. No RNA-only result here identifies that bridge in biological tissue.

Normal-only ridge fitting uses identities 0–39; identities 40–79 select from the
fixed 14-member path, including exact no-correction. Opposite binomial molecule
halves supply predictors and targets. With Gamma variation the halves share a
latent rate; marginal independence is explicitly not claimed. Both query halves
are predicted without query-dependent coefficient or penalty fitting. There is no
full-normal refit after selection.

## Locked experiment and independent tests

There are 32 independent base seeds, paired across seven observation conditions
and two depths: 448 datasets, **zero biological specimens**. Each contains 80 normal
and 80 query observations, 24 regions, 20 genes per region and four arms. Depth
1e6/1e7 is over the 480 retained genes, not whole-assay coverage. Half the query
population has a log(1.5) change on one arm, independently of state. This is a
high-count mechanism challenge, not a single-cell or whole-transcriptome realism
claim.

All observations precede fitting. The observation manifest has only case identity,
path and hash. A separate fitting process locks all 448 RNA predictions before a
separate scoring process reads any truth. Retained comparators include regional
count ratios, gene means/medians, leave-one-query-out profile/count pseudobulk,
the arm-high-pass ablation and stronger median/pseudobulk baselines using all
80 normals. Full-expression and truth-assisted oracle arms remain in the score
table; the oracle is not a usable competitor. All methods share feature support.

The 64 primary comparisons cross two state laws, two depths, two errors and eight
comparators. Each requires positive mean improvement and a simultaneous one-sided
exact 95% lower bound on strict-win probability above 1/2. Bonferroni protects this
64-comparison family. Rational binomial tails independently establish that 26/32
strict wins are required. Ties do not count as wins. Safeguards separately require
no point-mean error increase over raw gene means in clean, pure-dosage and
noise-only strata. These finite-simulation safeguards are not biological error
bounds or population noninferiority tests. Any failure rejects continuation under
the locked design.

The independent numerical and result checks pass:

- Noiseless dosage cancellation, state recovery, marker-free nonidentifiability
  and gene permutation fixtures pass. Augmented SciPy primal and kernel-dual ridge
  predictions differ by at most 6.67e-15 on the checked path points; condition
  numbers and operation/magnitude-based diagnostic allowances are retained.
- All 4,928 score rows retain 11 methods. A separate reader reconstructs all
  9,856 squared-error scores using extended-precision reductions; the maximum
  relative discrepancy is 3.39e-16. Every input, prediction and parameter hash
  verifies. All 279 selected no-correction outputs equal their raw profiles exactly.
- All 64 primary and 12 safeguard decisions agree with independent rational sums
  and integer binomial tails. Query reversal and count perturbation leave normal
  model selection, full selection-loss paths and normal identities unchanged in
  the two fixed isolation cases.

Numerical allowances concern arithmetic checks only. They do not relax any
scientific gate. The six failed mean differences range from 5.11e-9 to 2.61e-6,
far above the observed score-reconstruction differences. No numerical ambiguity
was used to convert a failed decision to a pass.

## Complete outcome and failure mechanism

All 64 primary comparisons have 32/32 wins; the smallest simultaneous lower
probability bound is 0.7996495, conditional on these generating laws. Selected
whole/cell-specific mean squared errors are:

| State law / retained-gene depth | Raw gene mean | Within-region correction | All-normal count pseudobulk |
|---|---:|---:|---:|
| Shared / 1e6 | 0.0503255 / 0.0483684 | 0.000080303 / 0.000077651 | 0.00872818 / 0.00790221 |
| Shared / 1e7 | 0.0502342 / 0.0482761 | 0.000007958 / 0.000007649 | 0.00872558 / 0.00790080 |
| Shared + Gamma noise / 1e6 | 0.0531098 / 0.0510720 | 0.00289968 / 0.00280014 | 0.00874399 / 0.00790732 |
| Shared + Gamma noise / 1e7 | 0.0528827 / 0.0507987 | 0.00279500 / 0.00269331 | 0.00873482 / 0.00790445 |

This mechanism removes marked, transferable state while preserving the planted
arm contrast around log(1.5). The full-expression ablation reduces that contrast
to about 0.319 in the Poisson state challenges. The comparison supports the
within-region restriction under this generator; it is not a proof for biological
CNA response, arbitrary gene mixtures or individual-cell DNA.

The six failed safeguards must remain explicit:

| Condition | Depth | Error endpoint | Correction minus raw mean MSE |
|---|---:|---|---:|
| Clean | 1e6 | Whole | +2.05722417e-8 |
| Clean | 1e6 | Cell-specific | +2.16207007e-8 |
| Clean | 1e7 | Whole | +5.46306350e-9 |
| Clean | 1e7 | Cell-specific | +5.11449954e-9 |
| Pure dosage | 1e6 | Cell-specific | +2.61395685e-6 |
| Pure dosage | 1e7 | Cell-specific | +2.07086742e-7 |

[Inference] Clean normal-validation noise sometimes selects a finite correction:
16/32 and 18/32 seeds at the two depths. In pure-dosage queries, the reference
contains real state and selects a nonzero predictor, but query sampling noise
still perturbs its inferred state. Whole-profile reference bias improves while
cell-specific error increases. Ideal mean-level protection is therefore weaker
than an observed-count guarantee. A query-specific abstention rule would be a new
method requiring a new frozen study; no such rule was tuned on these outcomes.

The reversed-transport challenge is consequential: at depth 1e6, whole MSE rises
from raw 0.0506954 to corrected 0.195009. Unmarked state remains unresolved, around
0.05056 versus raw 0.05033. The former falsifies unrestricted transport; the latter
is consistent with the independently demonstrated identification limit. More
optimization or count depth does not supply missing biological correspondence.

## Compute, provenance and decision

Both jobs use nslab account/QOS, CPU partition, Intel v4, one CPU, 8 GB and one-hour
limits, with existing Conda Python and isolated imports. They use no GPU, array,
interactive QOS, new package installation or modified shared environment.

| Job | Work | Completion | Elapsed | Peak process RSS |
|---|---|---|---|---:|
| 23205822 | Numerics, generation, fitting, scoring | 448/448 ██████████ | 2m13s | 93,876 KiB |
| 23205989 | Independent reader and isolation | 76/76 ██████████ | 11s | 62,668 KiB |

Running 0; pending 0; completed 2; failed jobs 0. The reader's afterok dependency
succeeded. Scheduler RSS sampling records 72,728 KiB for the main batch and zero
for the short reader; zero is not a measured zero-memory requirement. All four
main-process resource logs are retained separately. A preliminary login static
check used unavailable `python`; its AST step was not executed. Repeating with
`python3` passed, and every executable scientific check ran on compute nodes.

Runtime: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1; the main job ran on
ne1dc4-003. Source identities, environment/NumPy build configuration, normal roles,
full selection paths, input/prediction/truth hashes and exact submissions are
recorded. All 91 validated nonbinary package/data/test identities and all 276
preceding ledger entries remained unchanged before this entry. The existing
521-Python/four-Rust validation remains the relevant package check; it was not
rerun because production code and tests did not change.

**Decision:** retain the research implementation and conditional gain, reject its
general correction gate, and do not expand onto real data under this design.
Original production recovery remains 2/5; the continuous-state false-split result
remains 14/200; preserved correlated selected-interval coverage remains 54/1,000.
No biological calibration, uncertainty improvement or new clone recovery is claimed.
Independent corresponding DNA plus replicated RNA/reference transport evidence
would be required to extend the scientific interpretation. Biological replication
and suitable membership truth remain evidence gaps.

Evidence chain:
[frozen design](../dev/within_region_nuisance_20261110/DESIGN.md),
[estimator](../dev/within_region_nuisance_20261110/estimator.py),
[independent generator](../dev/within_region_nuisance_20261110/generator.py),
[executed workflow](../dev/within_region_nuisance_20261110/study.py),
[main submission](../dev/within_region_nuisance_20261110/submission.json),
[source identities](../dev/within_region_nuisance_20261110/source.sha256),
[full report](../dev/within_region_nuisance_20261110/run_23205822/results/REPORT.json),
[every score](../dev/within_region_nuisance_20261110/run_23205822/results/scores.tsv),
[all means](../dev/within_region_nuisance_20261110/run_23205822/results/means.tsv),
[numerical references](../dev/within_region_nuisance_20261110/run_23205822/results/NUMERICS.json),
[environment](../dev/within_region_nuisance_20261110/run_23205822/results/environment.json),
[independent reader](../dev/within_region_nuisance_20261110/read_results.py),
[reader submission](../dev/within_region_nuisance_20261110/reader_submission.json),
[all recounted decisions](../dev/within_region_nuisance_20261110/reader_23205989/RECOUNT.json),
[current source/ledger preservation](../dev/within_region_nuisance_20261110/current_source_check.json),
[all job outcomes](../dev/within_region_nuisance_20261110/accounting_complete.tsv),
[per-process resources](../dev/within_region_nuisance_20261110/resources_complete.json),
[final source integrity](../dev/within_region_nuisance_20261110/source_integrity_complete.json).
