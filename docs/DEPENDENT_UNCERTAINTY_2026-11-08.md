# Dependence-aware protection of selected RNA intervals

The objective is to narrow systematic-RNA intervals without weakening coverage
after adaptive interval selection. The new research calculation accounts for
overlap and covariance among all eligible intervals. It improves some simulated
power measurements, but does not resolve the missing biological covariance and
reference information. Production inference and uncertainty defaults are unchanged.

## Estimand, assumptions and construction

The estimand remains the interval mean systematic log-intensity effect under
Y_j given Z distributed Poisson(E_j exp(mu_j+Z_j)), with Z Gaussian and covariance
Sigma. E, Sigma and the bin/gap universe must be fixed independently of the
observations. The complete family contains1,860 intervals over120 bins in four
blocks. Regional transcription and reference mismatch remain alternatives to a
DNA alteration; the calculation cannot resolve that biological ambiguity.

Each independent calibration draw uses the maximum absolute standardized sum
over the entire interval family. This directly represents dependence among
overlapping intervals. It does not select an observed promising interval or fit
covariance to tumour residuals.

The total error budget remains0.05. Count bounds retain0.025. The remaining
budget is divided equally between a future Gaussian-tail event and calibration
failure,0.0125 each. This is a prespecified conservative allocation, not a
data-optimized operating point. Three covariance settings share the calibration
failure budget,0.0125/3 per bank.

With20,000 calibration maxima, choose the smallest one-based order-statistic rank
k satisfying P[Binomial(20000,0.9875)>=k]<=0.0125/3. This gives k=19791 and failure
probability0.00413578 per bank. The independent beta-order-statistic calculation
agrees. On successful calibration, the conditional future tail is at most0.0125;
including bank failures and count errors gives the stated model-conditional
repeated-calibration-and-data bound. It is not a deterministic guarantee for every
fixed bank. This distinguishes confidence in population content from coverage
itself, as in [NIST's tolerance-limit description](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/tolelimi.htm).

The final bank size is fixed before outcomes. Prefixes of5,000 and10,000 draws
are convergence diagnostics only. They do not choose the operating threshold.
Independent checks cover complete interval enumeration, gap boundaries,
translation/scaling behavior of standardized sums, positive covariance-factor
reconstruction, binomial/beta rank identity and a one-interval normal reference.
Floating-point sum error is padded conservatively. The theoretical guarantee
assumes the Gaussian calibration law; numerical checks are not a proof of exact
arithmetic or of the biological observation model.

## Frozen-panel results

All10,800 original intervals and their selections are retained. No HMM refit,
new simulated specimen or changed effect size enters the comparison. The six
candidate coverage endpoints join the six original simultaneous endpoints in a
12-member family, with the same95% target, exact undercoverage tests with Holm
adjustment, and simultaneous Clopper–Pearson intervals. Earlier legacy-SE failures
remain in their original reports.

All six candidate endpoints cover1,000/1,000 null simulations. Their simultaneous
coverage lower bound is0.993845; the original method has the same counts. These
are conditional simulation results, not biological replication.

| Condition | Original critical value | New critical value |
|---|---:|---:|
| Homogeneous independent bins | 4.35280 | 4.28689 |
| Heterogeneous independent bins | 4.35280 | 4.24001 |
| Heterogeneous correlated bins | 4.35280 | 3.96623 |

For the fixed six-bin interval at the original gain RNA effect0.243279,
zero-exclusions change from159 to161 of200 homogeneous simulations,3 to4 under
heterogeneity, and1 to3 under correlated heterogeneity. At twice that RNA effect,
the correlated condition changes from14/200 to30/200; its selected-interval
candidate has29/200. All power comparisons are secondary descriptive operating
measurements, not new confirmatory superiority tests. At the larger correlated
effect, median width changes from1.4054 to1.2921 natural-log units. There are no
infinite intervals in this positive-count panel.

**Decision:** the dependence-aware calculation provides a conditional improvement,
but power at the smaller effect remains poor. Keep it as a research reference.
It does not justify automatically replacing routine SEs or claiming confidence
for genomic CNA states. Before practical promotion, independently identified
reference and covariance inputs must exist, and repeated-workflow coverage must
include their estimation uncertainty. More simulations alone cannot supply them.

## Execution and evidence

Job22964546 completed the first run. Its raw log undercoverage probability for
complete coverage rounded slightly above zero. Job22964572 limits that report to
the exact probability range; the scientific intervals and calibration artifacts
are byte-identical between runs. Both outputs remain preserved. The second run
took12s with process peak RSS162MiB; audit22964574 completed in5s.

The audit verifies six recorded input hashes,10,800 unchanged original intervals,
10,800 unchanged candidate selections, all calibration ranks, scalar maxima for
nine prespecified draws, all12 coverage counts/bounds and48 power endpoints.
It also verifies395 preserved baseline files and65 unchanged package/test/native
files. The unchanged full package suite was not rerun.

Implementation: `benchmarks/hardening/gaussian_interval_maximum.py`,
`dependent_uncertainty_panel.py`, `dependent_uncertainty_audit.py`.
Immutable source snapshots, design, independent calibration banks, output files,
wheel identity and scheduler/resource records are under
`dev/dependent_uncertainty_20261108`.
