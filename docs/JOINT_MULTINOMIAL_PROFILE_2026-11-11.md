# Shared RNA populations:high-depth mechanism gain,inadequate general profile model

The scientific goal remains accurate RNA-derived genomic profiles and supported
CNA populations/spatial assignments,independently validated at the claimed
resolution. The new joint RNA model passes18/32 prespecified comparisons and
fails its overall gate. It captures a clean high-depth contrast but does not
solve low-depth whole-profile error or continuous state. Reject promotion and
retain all negative results. No existing package output or criterion changes.

## Assumption,estimand and exact model

Earlier independent regional priors cannot share a cell's repeated gain/loss
information across regions. The smallest joint alternative models each cell's
complete24-region count vector as one or two multinomials,conditional on its
observed total. Templates and mixture weights maximize a declared MAP objective:
summed log mixture count scores plus0.5 times all template/weight log probabilities.
These are Dirichlet parameters1.5,not an exact uncertain-membership posterior.
The component coefficients cancel in comparisons. Closed M steps add0.5 to
responsibility-weighted counts and component masses.

The estimand is regional relative RNA log odds. The action averages logit(template)
using the target's RNA component responsibilities. Two-component fitting uses
exactly collapsed and leading-Pearson-singular-vector starts,200 EM steps per
start,and the larger declared penalized objective. Every step checks monotonicity
within operation/magnitude roundoff. A400-step diagnostic cannot replace the
primary;wins must hold under BOTH budgets. Neither budget certifies convergence
or a global optimum. RNA components are not validated genomic clones.

All outer target cells are excluded from fitting AND model selection. Two internal
cell folds,defined by immutable ID parity,select one versus two components by
summed held-out RNA predictive score. This is not a genomic test,Bayes factor,
BIC or calibrated p-value. The always-two control isolates model-selection loss.
No DNA,truth label,larger K,new state factor,prior or favourable fold selects
the candidate. Same original counts,reference and two outer folds are retained.

Identification requires an adequate multinomial observation law and transport
of training RNA populations to targets. Continuous or discrete transcriptional
state,reference mismatch and biological sampling heterogeneity remain alternatives
to genomic dosage. Corresponding DNA and independent biological replication
must distinguish genomic explanations. This exposed generated panel has
biological n=0;48 base draws per condition are not patients.

## Computational evidence

All288 predictions lock before scoring truth. All1152 original baseline score
rows match exactly. Qualification includes rational one-component templates and
M steps,an exact255/257 two-component action,complete SciPy multinomial PMFs,
label invariance,exact ID permutation and outer-target perturbation exclusion
from templates and K selection. Independent SciPy optimization agrees on the
small fixture within1.77636e-15 against2.47050e-12 budget. That fixture does not
certify a global optimum for every panel fit.

Independent reader23210952 reconstructs3317760 action values using complete
multinomial PMFs,with maximum discrepancy1.01208e-12. It reconstructs counts,
2304 inner-selection checks,independent endpoint sums and all32 exact binomial
decision tails. Every original50%-win lower-bound criterion with Bonferroni32
remains required;all48 draws remain in each denominator. Both scientific
decision and failed baseline comparisons reproduce.

## Performance and prespecified ablations

Mean whole-profile RNA log-odds MSE;lower is better. Every row uses48 original
generated specimens and the same scoring support.

| Condition | Reads/cell | Pseudobulk | Moment-Beta expected logit | Joint primary | Always two control |
|---|---:|---:|---:|---:|---:|
| Two CNA populations | 300 | .004133 | .005170 | .004862 | .006244 |
| Two CNA populations | 3000 | .003485 | .001268 | .000330 | .000330 |
| Two populations plus state | 300 | .016197 | .011178 | .013154 | .012216 |
| Two populations plus state | 3000 | .015545 | .002314 | .008120 | .008120 |

[Observed pattern,not a guarantee] The clean high-depth joint model improves
both whole/minority error against every retained baseline and retains100.73%
of the generated minority contrast. This conditional gain is insufficient for
the overall32/32 gate. The clean low-depth primary wins0/48 whole-profile
comparisons against pseudobulk. Its always-two control ALSO wins0/48 and has
worse mean error;forcing a split is not a sufficient repair.

Training-only prediction selects two components in2/96 outer fits at low-depth
clean two-population counts,versus96/96 at high depth. At low-depth plus state,
it selects76/96;at high-depth plus state96/96. Homogeneous controls select one
inall96 fits at both depths. Counts are fit decisions,not biological n or
false-split estimates. Selection choices agree under both iteration budgets.

The clean low-depth always-two control improves minority mean MSE but still
fails the whole-profile criterion. At high depth with state,the selected model
retains only14.46% of the genomic contrast,and minority RNA MSE is.081232 versus
.014072 for the marginal Beta action. [Inference] Shared discrete templates
alone are inadequate for these continuous-state observations. These data do
not justify increasing K or adding a latent factor after the failed gate.

The independent ablation summary retains every method/condition mean and all32
always-two comparisons. It cannot promote the control in place of the primary.
No primary point-win sign changes under doubled iterations. This is a model/
information failure on these conditions,not a numerical decision ambiguity.
State-free genomic-response errors are retained separately from RNA errors.
Improving prediction of a transcriptional state would not validate genomic dosage.

## All execution and preservation

All three jobs finish0:0 on nslab/nslab CPU,AMD-v6,ne1dc6-009,using the original
read-only scientific-review-response Conda environment/runtime003 and one thread.
No new dependency,GPU,interactive QOS,array,shared installation or running-job
change. The main source_v1 was prepared but unexecuted;source_v2 includes a
pre-submission Python-bool serialization correction. The executed main manifest
is unchanged by the later separately frozen ablation summarizer.

| Job | Work | Requested CPU/RAM/time | Actual elapsed | Peak process RSS,KiB |
|---|---|---|---|---:|
|23210951|Qualification,288 fits,scoring|1/8GB/6h|2m55.15s|81756|
|23210952|Independent PMF/decision reader|1/4GB/1h|6.82s|82176|
|23211037|Prespecified ablation summary|1/4GB/1h|2.72s|12288|

[Frozen contract](../dev/joint_multinomial_profile_20261111/DESIGN.md),
[executed estimator](../dev/joint_multinomial_profile_20261111/model.py),
[qualification](../dev/joint_multinomial_profile_20261111/run_23210951/results/QUALIFIED.json),
[all32 primary comparisons](../dev/joint_multinomial_profile_20261111/run_23210951/results/REPORT.json),
[independent reader](../dev/joint_multinomial_profile_20261111/reader_23210952/results/REPORT.json),
[all ablations/means/selection counts](../dev/joint_multinomial_profile_20261111/ablation_23211037/results/REPORT.json),
[all scores](../dev/joint_multinomial_profile_20261111/run_23210951/results/scores.tsv),
[exact submission commands](../dev/joint_multinomial_profile_20261111/submissions.json).

Original production recovery2/5,state false splits14/200,selected coverage54/1000,
the5% simultaneous criterion and frozen artifacts remain unchanged. No biological
calibration,competitor superiority,default promotion or SOTA completion follows.
