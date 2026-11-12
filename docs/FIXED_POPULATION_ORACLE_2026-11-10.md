# Ideal-observation power does not establish an explanation for the recovery failures

Job23165792 completed0:0 in32s. The numerical checks pass. The simultaneous
95% Monte Carlo upper bound on the optimal5%-FPR ideal-observation power is
0.9913789356. This high ceiling does not establish a restrictive information
limit explaining the failed recoveries. It also does not prove that RNA counts
can attain that power or that the existing method is inefficient. Do not use
this calculation to excuse a regression or relax a scientific criterion.

The hypotheses are deliberately specific. Both use600 cells and the original
shift delta=.2432790649. H1 has exactly240 shifted cells and360 unshifted cells,
with independent Gaussian regional expression state of SD delta*sqrt(.24).
H0 has one CNA population and a Gaussian regional state whose mean and variance
match H1's marginal mean/variance. This is an additional plausible RNA-state
competitor, not a new interpretation of the original phase8 state-null errors.

The observer sees each cell's ideal regional log-rate before gene/count
sampling. Conditional on that rate and the common baseline/exposure, the
specified generator has the same RNA law in both hypotheses. Its RNA counts
therefore cannot contain more discrimination information than the ideal rate
for this pair of hypotheses. The bound applies to these fixed known parameters;
it is not a universal limit for other effects, observations or biological states.
In particular it does not validate genomic labels or same-cell membership.

The alternative likelihood averages over all choose(600,240) assignments. A
log-space elementary-symmetric-polynomial recurrence evaluates that fixed-count
law without enumerating every vector. An iid Bernoulli mixture would be a
different law and was not substituted. Small exhaustive label enumeration,
permutations, label complement, pure-label cases and eight full600-cell
extended-precision comparisons pass. The largest full600 numerical discrepancy
is6.2528e-13, within its magnitude/operation allowance4.1366e-10. This allowance
is checked arithmetic diagnostics, not directed-rounding certification.

There are10000 independent null calibration draws and20000 alternative draws,
with the frozen seed and128-row batching. The97.5% order-statistic bracket for
the null95th percentile uses ranks9451/9549 and gives critical scores
[-1.4069793920,-.9806714684], widened by arithmetic allowances. At the liberal
lower score,19801/20000 alternative statistics possibly exceed the threshold.
A one-sided exact97.5% binomial upper bound, combined with the quantile bracket
by the union bound, gives the recorded95% Monte Carlo upper power envelope.
Beta inversion and an independent binomial-CDF root agree exactly. The19801
count is not achieved pipeline power at5% false-split error.

All Monte Carlo populations are numerical replications; biological n=0. No
package inference, thresholds, frozen outputs, independent validation status
or promotion decision changes. The original failed candidate still fails its
false-split/power gate. The next useful work remains the strict refinement
intervention and exact tracing of the other under-split/over-split cases.

Process runtime26.25s, peak138304KB. The CPU/innovation allocation requested
4CPUs/8GB/1h; the validated environment and911-file snapshot were read only.
All909 prior source files remain exact. No new dependency or shared install.

[Design fixed before execution](../dev/fixed_population_oracle_20261110/DESIGN.md),
[executed calculation](../benchmarks/hardening/fixed_population_oracle.py),
[complete report](../dev/fixed_population_oracle_20261110/oracle_v1/report.json),
[complete manifest](../dev/fixed_population_oracle_20261110/oracle_v1/run_manifest.json),
[resource record](../dev/fixed_population_oracle_20261110/analysis_23165792/resources.txt),
[submission command](../dev/fixed_population_oracle_20261110/jobs.json).
