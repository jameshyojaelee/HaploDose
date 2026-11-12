# Observable-mixture confirmation and corrected profile units

The objective remains accurate RNA-derived dosage profiles and CNA-defined clone
recovery, ultimately tested against fixed strong comparator families on suitable
independent evidence. The preceding goal turn made progress through factorial
attribution and a stronger profile baseline. This turn tests an actual model
intervention and corrects a newly discovered baseline-comparison defect.
Performance findings below are [Observed pattern, not a guarantee].

## Conditional mixture confirmation does not repair the remaining failures

The candidate retains the existing event proposals and gene veto, profile
initialization, depth matching, posterior cutoff .95, minimum segment length,
and every original recovery assertion. It compares a latent log-rate state delta
with the expected contrast in mixed assigned pools:

`log(1+p_a*expm1(delta)) - log(1+p_b*expm1(delta))`.

Pool purity p is the proposal component probability averaged with the actual
pooling weights and cell exposures. Its allocation factor is
`h=sum((w*L)^2*q*(1-q))/sum(w*L)^2`. The candidate adds conditional delta-method
allocation variance `h*expm1(delta)^2/(1+p*expm1(delta))^2` from each disjoint pool.
Pure memberships reproduce the original state effects exactly; equal purities
require abstention. Both initial and depth-matched confirmation use the changed
observation model. The existing experimental merge weighting is held fixed in
both candidate and control.

The estimand remains an assigned-pool RNA contrast under a latent mixture working
model. Common bin composition after nuisance correction and accurate fitted
probabilities are identifying assumptions. Allocation variance does not include
probability/nuisance fitting uncertainty. Shared cell membership induces
cross-bin dependence not represented by the diagonal HMM emissions. Correlated
expression state can reproduce genomic-looking structure. No independent
confirmation, exact posterior calibration or biological validity is claimed.

Job 22961219 verifies means and variances against exhaustive enumeration of all
16 latent allocations for three state effects. Maximum moment discrepancy is
3.89e-16. Pure-membership identity and weight-rescaling invariance pass. These
checks validate the implemented conditional calculation, not its applicability.

Jobs 22961220 and 22961221 replay all five unchanged frozen cases. Both paths pass
20261031, 20261107 and 20261108 and fail 20261109 and 20261110: **3/5 versus 3/5**.
The remaining chr12 proposals still have no supported initial confirmation
segment. Their fitted low-component pool purities are approximately .684/.257
and .777/.228. Changing the expected mixed-pool state and adding allocation
variance is insufficient to pass the original recovery requirements.

This candidate fails the prespecified recovery-improvement gate. It is retained
as a negative experiment and is not promoted. The larger repeated-workflow
null/power sweep was not run: there is no demonstrated recovery benefit to
advance. The prior 6.81% false-split upper bound still fails the 5% criterion.
A next candidate must retain per-cell likelihood information and shared latent
membership across bins, rather than only repairing pooled means. Numerical
correctness and full-workflow calibration will remain separate requirements.

## Correction: TCGA shared segmentation penalty had incompatible units

The earlier operator experiment in `dev/attribution_20261108/profile_v1` treated
the TCGA naive track as natural-log RNA effects. Its source,
`benchmarks/tcga_bulk/subarm.py::type_windows`, explicitly constructs
`log2(CPM+1)`, while HaploDose's bulk log ratios use natural logs. Rank correlation
alone is invariant to the scale conversion; a shared quadratic-loss penalty is
not. The first TCGA shared-penalty comparison therefore did not match operators
in effect. **Its conclusion that the naive track wins that comparison is
withdrawn.** The original output and decision remain preserved with this erratum.

`profile_operator_control.py` now converts fitting values into natural logs and
records input and output units explicitly. Original input columns remain intact.
A test verifies that multiplying data by log(2) and the quadratic penalty by
log(2)^2 preserves partition boundaries and scales fitted values consistently.
No fitted parameter, scientific threshold or case selection changes.

Corrected job 22961363 gives TCGA Spearman agreement **.8038 naive versus .8074
HaploDose**, rather than the earlier .8038 versus .7847. This single transformed
input remains a processing diagnostic, not raw-count calibration or evidence of
population superiority. Both Kinker controls and both U-CAN results are unchanged:
naive is higher in the two Kinker lines, while HaploDose is higher in the U-CAN
tumours. No general superiority claim is supported by these selected cases.

## Preservation and execution

The dirty baseline and Git patches are captured under
`dev/mixture_confirmation_20261108/baseline`; every execution uses an immutable
source snapshot and the previously verified wheel read-only. All numerical
checks and replays run on Intel v4 compute nodes under innovation or nslab.
Production package code is unchanged; a new full package suite was not run.
The latest package suite remains 319 Python passes, one GPU-only skip and four
Rust passes. New source affects only experimental benchmark paths and reports.

Independent audit job 22961634 confirms both 362-file frozen source inventories,
unchanged original criteria/generators/package files, exact preservation of the
prior merge control, and **exactly identical final posterior arrays between
candidate and control in all five cases**. It also verifies that the four
already-consistent profile cases are unchanged and every original profile input
column is retained. All five jobs complete successfully; the peak scheduler RSS
is approximately 1.04 GiB. Short-job RSS reported as zero is missing accounting,
not zero memory use; per-job `resources.txt` also records process measurements.
