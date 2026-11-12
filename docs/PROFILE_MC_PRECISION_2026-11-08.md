# Numerical uncertainty in the existing profile bootstrap

Update: the fixed one-million-draw refinement and independent audit now resolve
the correlation calculations' numerical precision. The original20000-draw
finding below remains as history. No biological/reference limitation is removed.

The original62 correlation results reproduce exactly before extension. All93
scientific endpoints remain;31 loss endpoints retain their original precision.
For the programme-direct minus programme-median Pearson contrast:

| Resampling draws | Holm p | Basic simultaneous interval |
|---:|---:|---:|
| 20000 | .0115494 | [.00737175,.06252234] |
| 200000 | .0192499 | [.00695659,.06283618] |
| 1000000 | .0186340 | [.00644310,.06280335] |

At the fixed final budget its simultaneous numerical Holm upper bound is.0234591.
The numerical band for the lower interval endpoint is[.00568933,.00700378];
for the upper endpoint it is[.06239981,.06331007]. All thirteen originally positive
correlation claims satisfy the numerical stability checks at the final budget.
This supports the narrower numerical interpretation of those conditional
development comparisons, not independent superiority or certified bootstrap
coverage. The new method still does not establish a benefit over the original
simple baseline or a benefit from adding programme transfer itself.

Refinement22997407 completes in1m47s; independent audit22997510 in46s. Process
peak RSS is640328KiB and635588KiB respectively. The audit independently checks434
scalar draw endpoints,62 probability counts and248 interval endpoints at each
checkpoint, the full Holm family, unchanged scientific populations and unchanged
loss results. It uses exact rational interpolation of binary64 order statistics
with an operation-derived rounding allowance. All million draws remain saved.
Only the final fixed budget carries the stated numerical confidence; earlier
checkpoints are convergence diagnostics. A syntax-error attempt22997381 stopped
before analysis and remains preserved. Sources refinement_source_v2 and
refinement_audit_source_v1 produced refinement_v2 and refinement_audit_v1.

## Original finite-draw finding

The20000-draw profile report has insufficient Monte Carlo precision to certify
all of its extreme interval endpoints at the additional simultaneous numerical
confidence assessed here. One programme-transfer comparison also has an adjusted
probability uncertainty range crossing0.05. This is a numerical precision gap,
not evidence that its biological effect is zero. Original results are preserved.

The calculation conditions on the same patients, thresholds, profiles and finite
eight-reference-refit ensemble. It reconstructs integer exceedance counts from
the reported plus-one probabilities exactly. A separate95% numerical family
contains93 exceedance probabilities and186 interval quantiles,279 in total.
Exact binomial bounds are propagated through monotone Holm adjustment. These
bounds describe resampling error, not patient uncertainty or bootstrap validity.

At the scientific interval tail probability0.05/(2*93), only5.376 bootstrap draws
are expected. The simultaneous numerical order-statistic ranges are ranks0–17
and19984–20001. Rank0 means negative infinity and rank20001 positive infinity;
neither is replaced by a convenient finite extrapolation. Thus those numerical
endpoint bands are unbounded at this diagnostic confidence level.

For Pearson programme_gene_partition minus programme_median_partition, the
original exceedance count is2/20000 and reported Holm p is0.0115494. Its
simultaneous Monte Carlo Holm upper bound is0.0540944. The numerical precision
of that decision is unresolved. The other twelve original superior flags have
Monte Carlo Holm upper bounds below0.05, but this does not resolve their interval
endpoint precision or validate biological coverage.

The next admissible correction is a fixed-budget numerical refinement of the
same bootstrap law, with original-draw reproduction and convergence checks.
It must not change methods, comparator families, thresholds or biological n, or
sample until a preferred conclusion is obtained. No inference default is changed.
Independent reference/covariance knowledge and DNA truth remain separate gaps.

Job22997165 completes in9s with process peak RSS119844KiB. Source candidate_v2,
output audit_v2, is under `dev/profile_mc_precision_20261108`. Job22997145 and
candidate_v1 are preserved: that attempt expected a manifest the older audit
did not write. The corrected reader verifies its actual zero exit, executed
script/source identity and equality of the report with captured job output.
All93 comparisons remain. No additional inference or biological data was used.
