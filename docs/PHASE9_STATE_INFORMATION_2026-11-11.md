# Phase9 state and population information:what the oracle can distinguish

The specified phase9 amplitudes permit high population detection power when
observed without count noise and with known distribution parameters. At a fixed
5% specimen false-positive rate,the optimal likelihood-ratio test has power
bounded by **0.9861029190–0.9911421822** with at least95% joint statistical
confidence,conditional on the checked numerical enclosures. This does not
establish practical RNA power,genomic identity,or error control of HaploDose.

The independently checked result rules out importing the earlier low-information
balanced experiment's bound into this different generator. It leaves count
noise,nuisance estimation,conditional scoring and fitting as explanations for
the actual candidate's lost true contrast. No frozen result or criterion changes.

An earlier[fixed-population calculation](FIXED_POPULATION_ORACLE_2026-11-10.md)
already reported a0.9913789356 upper bound for the same underlying600/240
effect/state law in phase8. The new computation is an independent replication
with an added lower confidence bound,exact rank-tail verification and stronger
power-inversion checks. It is not a new discovery of high available information.
This precedent was missed when preparing the new run. No further latent-oracle
replication is warranted;the next experiment uses the saved count measurements.

## Estimand,observation model and identifying restrictions

The estimand is the optimal specimen-level rejection probability under two
specified simple probability laws,at type-I error0.05. Observations are600
noiseless regional RNA log-amplitudes x,not counts or DNA. Alternative P has
exactly240 uniformly permuted high cells,within-class Gaussian SD
sigma=0.11918191480409052,and separation delta=0.24327906486489861. Withw=.4,
component means are−w*delta and(1−w)*delta. The null Q has one Gaussian RNA
population with mean0 and variance sigma²+w(1−w)delta². It matches the alternative's
marginal mean/variance but not its full distribution or fixed-count dependence.

A fixed regional offsetw*delta is shared under both hypotheses. The original
count kernel would receive eta=x+w*delta in either case. This is a common
coordinate change,not hypothesis-specific fitted centering. Q is a broader state
alternative than the original phase9 zero-mean/fixed-variance null. The original
14/200 state-error result and its11.87147% simultaneous upper bound remain.

The frozen generator and an archived generation record verify n,k,delta,sigma.
Library exposures and baseline rates are generated independently of clone labels;
the count kernel depends on those labels through total regional amplitude.
Under these restrictions,counts and common independent nuisance measurements
cannot improve on the known-parameter amplitude oracle. High oracle power is
therefore compatible with low count-based power;it does not demonstrate an
available estimator or locate the practical information loss by itself.

Biological alternative:two transcriptional states could produce the exact same
amplitude mixture as two CNA populations. This diagnostic cannot distinguish
those mechanisms. Corresponding independent DNA,not the likelihood-ratio result,
would establish genomic populations. Biological n=0. The20000 specimens per
condition are independent simulations under the declared laws.

## Exact density and prespecified error/power calculation

The P density averages Gaussian products over every fixed240-of600 assignment.
A log-domain elementary-symmetric-polynomial recurrence evaluates the required
degree240 coefficient without choosing labels or approximating them as iid.
The log ratio to Q is the Neyman–Pearson statistic for these known simple laws.
There is no fitted HaploDose method,optimizer,region search or tuning step.

Two fixed seed sequences generate exactly20000 Q and20000 P specimens. The
95th Q quantile is bracketed by order statistics18930 and19069. Exact integer
binomial CDF recurrences verify their respective tail probabilities
0.01168302444 and0.01247218772,each below1/80,and verify that moving either rank
inward would fail that allocation. Outward numerical rounding gives threshold
limits−1.32838506816 and−1.08429151993.

There are19759 guaranteed and19792 possible alternative rejections at the
unknown true Q quantile. One-sided Clopper–Pearson lower/upper limits,each with
failure allocation1/80,produce the quoted power bracket. Together the four
failure allocations give at least95% joint confidence. Power-bound inversion is
independently enclosed with80- and100-digit binomial CDF recurrences and outward
binary64 rounding. No universal power threshold is introduced.

The numerical allowance counts recurrence operations and magnitudes;it is a
checked estimate,not formal interval arithmetic. Statistical confidence is
conditional on these numerical enclosures. Exhaustive unbalanced small-n SciPy
densities,zero-effect identity,600-cell permutation/extended-precision checks,
and exact agreement with the balanced implementation all pass. Eight full-size
independent posterior-count polynomial convolutions cover extremes and
threshold-near observations. Their discrepancies range below5.0e-13,against
operation/magnitude allowances around6.7e-9–7.4e-9. Every saved amplitude,label,
score and metadata hash verifies;all alternatives have exactly240 high labels.

## Compute,evidence and next decision

Jobs23211374/23211375 complete0:0 in1m18s/1m15s,Slurm peak88860/88856KiB.
Independent reader23211376 completes0:0 in6.17s,time-v peak236628KiB.
All use nslab/nslab CPU Intel-v4,1CPU4GB1h,existing isolated Conda environment,
neutral directories,one thread and before/after source hashes. No production
source,shared installation,frozen result or scientific criterion changes.

[Frozen design](../dev/phase9_state_information_20261111/DESIGN.md),
[density implementation](../dev/phase9_state_information_20261111/model.py),
[independent reader](../dev/phase9_state_information_20261111/reader.py),
[qualification checks](../dev/phase9_state_information_20261111/null_v1/CHECKS.json),
[complete independent result](../dev/phase9_state_information_20261111/reader_23211376/results/REPORT.json),
[all exact submission commands](../dev/phase9_state_information_20261111/SUBMISSIONS.json).

The next smallest experiment uses the SAME exposed paired state/true count
datasets and original event genes,with simulator-known count nuisance,unchanged
fit/score/numerical rules and minimum30. It tests whether nuisance specification
alone can rescue the paired distinction. This is an explicitly oracle diagnostic;
success would require subsequent estimation-focused ablations. Failure would
falsify this proposed sufficient correction,without proving RNA contains no
information. Jobs23211662/23211663 and independent23211664 are submitted on
nslab/nslab CPU Intel-v4,1CPU8GB3h/3h/1h;results pending.
[Exact follow-up contract](../dev/state_nuisance_oracle_20261111/DESIGN.md),
[follow-up submission commands](../dev/state_nuisance_oracle_20261111/SUBMISSIONS.json).

Full recovery,complete-workflow false-split/power evaluation,uncertainty
calibration and independent biological validation remain separate requirements.
