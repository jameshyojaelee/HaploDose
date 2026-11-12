# Why conservative RNA-effect intervals have low power

The goal is useful uncertainty for the systematic regional RNA effect, while
retaining coverage after interval selection. This analysis separates counting
noise, latent/reference variation and the simultaneous penalty. It creates no
new interval method, cutoff or biological validation claim.

The frozen model is Y_j given Z distributed Poisson(E_j exp(mu_j+Z_j)), with
independently fixed E and known Gaussian covariance Sigma. The target is the
mean systematic mu over the reported interval. Its alternatives include regional
transcription and reference mismatch; RNA intervals alone cannot identify DNA
dosage. Routine fitted residual covariance does not satisfy the known-Sigma
contract merely because the software accepts a matrix.

Job22964032 reconstructs all10,800 saved simultaneous intervals from the earlier
coverage and power experiments. No new draws, HMM fits, selected intervals or
scientific criteria are introduced. Independent scalar sums reconstruct their
widths to within6.66e-16. Normal-tail calculations are independently checked by
numerical quadrature.

| Fixed six-bin interval, null simulations | Median count contribution | Latent/reference contribution | Latent fraction of total width |
|---|---:|---:|---:|
| Homogeneous independent bins | 0.1665 | 0.2843 | 63.1% |
| Heterogeneous independent bins | 0.1668 | 0.6914 | 80.6% |
| Heterogeneous correlated bins | 0.1661 | 1.2752 | 88.5% |

Widths are in natural-log RNA-effect units. Each row uses the original1,000
simulations, not biological specimens. The current Gaussian critical value is
4.3528 because the procedure protects all1,860 eligible intervals and retains
its original split of the error budget between count and latent uncertainty.

At the original gain-state RNA separation0.243279, removing all counting noise
would leave the following fixed-interval exclusion probabilities under the model:

| Condition | Current all-interval latent bound, infinite depth | Known-covariance fixed-interval reference |
|---|---:|---:|
| Homogeneous independent bins | 99.90% | approximately100% |
| Heterogeneous independent bins | 9.86% | 86.50% |
| Heterogeneous correlated bins | 0.36% | 38.26% |

The second column is an infinite-depth limit, not a statement of monotone power
at every finite count realization. The final column protects one prespecified
interval only; it is invalid as a replacement for selection-aware inference.
It shows the cost of simultaneous protection under these assumptions. Even that
oracle reference has limited power in the correlated condition at this effect.

**Decision:** deeper counts alone cannot resolve the correlated-condition failure.
A dependence-aware simultaneous calculation is a justified bounded follow-up,
but any such improvement must retain the complete interval family, independently
verify numerical error and satisfy the same coverage/power experiment. It would
still require independently identified covariance and expectation for biological
use. More simulations cannot provide those measurements. The current routine SE
remains an uncalibrated scale diagnostic; the optional conservative helper retains
its explicit model contract and original implementation.

Implementation is `benchmarks/hardening/uncertainty_width_attribution.py`.
Immutable source, design, input/output hashes, execution environment and resource
records are under `dev/uncertainty_attribution_20261108`. The CPU job completed in
22s with process peak RSS149MiB. Frozen coverage failures and operating curves
remain unchanged; production, original tests and scientific thresholds are untouched.
