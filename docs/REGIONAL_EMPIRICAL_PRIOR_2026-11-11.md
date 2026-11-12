# Regional empirical prior: qualified computation,failed performance gate

The scientific goal remains accurate RNA-derived genomic profiles and supported
clone/spatial inference at the claimed resolution. This experiment changes the
regional prior distribution rather than repeating rejected posterior-action,
trunk-size or covariance changes. It does not pass its frozen scientific gate.
No package default or historical result changes.

## Estimand,model and fixed decision

The estimand is individual-cell regional RNA log odds relative to the specified
reference. A binomial RNA margin is mixed over a finite empirical training grid.
For each region and target-excluded fold,the support is unique
`(y+0.5)/(N+1)` training fractions plus the pooled training mean. The weight prior
is Dirichlet`1+1/M`:one total MAP pseudocount,not concentration one. Equality-
constrained damped Newton solves the strictly convex weight objective. Primary
action is posterior expected logit. A matched moment-Beta expected-logit control
separates a distributional change from an action change.

Identifying assumptions are exchangeable training/target probabilities and an
adequate binomial margin. Regional transcriptional state can produce the same
RNA profile as genomic dosage. State-free genomic-response error is separately
reported;better RNA prediction alone would not validate DNA or membership.
Independent corresponding DNA and biological replication are still missing.

The existing288 generated conditions contain48 base draws per scenario,paired
depth300/3000,120cells and24regions. These are exposed regression simulations,
biological n=0. The entire original1152 baseline score rows must match exactly.
The primary32 comparisons retain the earlier24 and add the matched-action
control:two heterogeneous scenarios,two depths,two endpoints and four baselines.
Every one-sided exact lower win-probability bound,Bonferroni across32,must exceed
0.5. All48 draws stay in every denominator;unresolved numerical boundaries do
not count as wins. No support,strength,action,fold or depth sweep after failure.

## Executed result

Study23210503 and independent reader23210504 complete0:0 on ne1dc6-009.
All1152 old rows match exactly. All288 predictions lock before truth scoring.
Qualification covers exact rational posterior identities,an asymmetric Fraction
root,independent SciPy optimization,finite-difference truncation/roundoff bounds,
permutations and held-out perturbation. The reader independently reconstructs
all829440 posterior actions from complete binomial PMFs;maximum discrepancy is
6.21725e-15. It checks inputs/folds/support and recounts all32 decisions with
exact rational binomial tails. Computational qualification passes.

The candidate passes27/32 comparisons. These five fail:

| Scenario | Depth | Endpoint | Comparator | Resolved wins/48 | Mean comparator MSE minus candidate MSE |
|---|---:|---|---|---:|---:|
| Two CNA populations | 300 | Whole RNA | Pseudobulk | 0 | -0.00325721 |
| Two CNA populations | 300 | Whole RNA | Global shrinkage | 0 | -0.00263022 |
| Two CNA populations | 300 | Whole RNA | Moment-Beta expected logit | 0 | -0.00222007 |
| Two populations plus state | 300 | Whole RNA | Moment-Beta expected logit | 0 | -0.00229097 |
| Two populations plus state | 3000 | Whole RNA | Moment-Beta expected logit | 25 | +0.00000104 |

Minority-region comparisons pass,but that cannot replace the failed whole-profile
criteria. Point wins and numerically resolved wins agree;this failure is not an
optimizer rounding ambiguity. Reject the unchanged prior. Variable-depth or
biological transfer of this failed candidate is not authorized by its gate.

## Compute,commands and artifacts

Both jobs use nslab/nslab CPU,AMD-v6,the original read-only scientific-review-
response Conda environment/runtime003 and one numerical thread. Study requests
1CPU8GB6h,finishes25.81s,peak81764KiB. Reader requests1CPU4GB1h,finishes6.54s,
peak82944KiB. Actual NumPy2.4.6/Python3.11.18 and node are recorded.

The exact sbatch argument vectors,UTC submission times and dependencies are in
[submissions.json](../dev/regional_empirical_prior_20261111/submissions.json).
The frozen [design](../dev/regional_empirical_prior_20261111/DESIGN.md),
[source identities](../dev/regional_empirical_prior_20261111/source_manifest.json),
[study report](../dev/regional_empirical_prior_20261111/run_23210503/results/REPORT.json),
[independent report](../dev/regional_empirical_prior_20261111/reader_23210504/results/REPORT.json),
[all scores](../dev/regional_empirical_prior_20261111/run_23210503/results/scores.tsv),
[study resources](../dev/regional_empirical_prior_20261111/run_23210503/resources.txt)
and [reader resources](../dev/regional_empirical_prior_20261111/reader_23210504/resources.txt)
preserve the complete negative result. No production promotion,independent DNA
validation,5% biological error guarantee or SOTA claim follows.
