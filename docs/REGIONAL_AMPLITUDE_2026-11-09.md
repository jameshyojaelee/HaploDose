# Shared regional RNA effects increase support but do not rescue confirmation

The fixed-state family was a remaining explanation for rejecting informative
proposals after pool arithmetic and nuisance removal passed their checks. This
experiment integrates one continuous RNA effect across each nonzero region.
It increases conditional support for the two rejected production contrasts,
but neither reaches the unchanged posterior/core requirements. Do not promote
this model or run a large null/power study as though recovery improved.

## Scientific contract and prespecified falsifier

The estimand is an assigned-pool regional RNA log-rate contrast, conditional on
the executed membership, nuisance correction, variance and count split. It is
not the difference between pure genomic clones. The observation model remains
Gaussian for the saved pooled contrasts. Within each maximal nonzero state run,
one shared effect has a normal prior centred at plus/minus the existing
beta*log(1.5) step, with SD equal to that step and truncated to its declared
sign. The zero state remains zero. This prior scale was fixed from the existing
effect scale before outcomes; no grid, truth fitting or selected width is used.

Conditional diagonal sampling variation and a constant regional effect are
identifying working assumptions. Misallocated memberships, correlated expression
state, reference mismatch and dependence remain alternatives. A shared RNA
effect can arise without a DNA alteration. Independent corresponding DNA and
biological replication would distinguish genomic meaning; numerical integration
cannot. Biological n=0. Both bin and known-size gene paths, all five exposed
canonical cases and all 27 saved binary comparisons remain in the record.

Under the proposed amplitude-restriction explanation, integrating a coherent
weak effect should raise regional support and might restore a qualifying core.
Failure to restore a core rejects this fixed candidate as a sufficient rescue
of those comparisons. A restored initial core alone would still not establish
final clone recovery, merge reliability or false-split control. Retain .95
support at every bin in a credible core and each comparison's original minimum
core length. No threshold or original recovery assertion changes.

The [original design](../dev/regional_amplitude_20261109/DESIGN.md) precedes all
new runs. Earlier [mixture-mean confirmation](MIXTURE_CONFIRMATION_2026-11-08.md),
[weighted stages](GENE_THINNING_WEIGHTED_2026-11-09.md) and
[executed-model reconstruction](CONFIRMATION_OBSERVATION_2026-11-09.md) remain
separate experiments. This does not relabel them as successes or repeat a
finer independent-bin effect grid.

## Exact model and independent checks

[Model implementation](../benchmarks/hardening/regional_amplitude_model.py)
integrates the regional scalar effect analytically. It sums and maximizes over
all state paths using a semi-Markov representation: diagonal transition
probabilities supply run lengths, and off-diagonal probabilities supply state
changes. The original state-path prior and genomic resets remain unchanged.
The integration is joint across bins in a run, not independent per-bin draws.
The model has quadratic block complexity; no new production dependency is used.

[Probe and audit](../benchmarks/hardening/regional_amplitude_probe.py) establish:

- Two exhaustive 243-path fixtures, including resets, agree with independent
  scalar integration. Posterior discrepancies are 8.88e-16 for the point prior
  and 5.55e-16 for the continuous prior. Corresponding checked guards are
  1.21e-12 and 2.76e-10. MAP paths, effect-sign reversal and common unit scaling
  agree. Scalar integration requests at 1e-9 and 1e-12 must converge under
  reported error plus accumulation estimates, rather than assuming a requested
  tolerance establishes accuracy.
- All 27 point-prior limits agree with independent long-double ordinary-HMM
  calculations. This checks equivalent mathematical models; it is not a claim
  of exact float32 native-array equality. The original native arrays were
  independently checked in the preceding observation study.
- All 27 candidate credible-core tables are independently recounted from saved
  posterior/path arrays under the original rules.
- The prespecified production109 first-half comparison is checked on its
  complete 19-position chr12 block, using independent amplitude quadrature
  represented as an expanded ordinary HMM. A microscopic amplitude persists
  on self-transition; a new run draws a new amplitude. Quadrature weights are
  not normalized or clipped. Gaussian prior tails, the maximum likelihood ratio
  and the all-null path bound the omitted mass.

The full-block quadrature changes across orders 32/64/128/256/512 are,
successively, 0.02210, 3.57e-5, 1.50e-11 and 2.23e-14. The final reference
differs from the analytical model by 2.87e-14, within the operation/magnitude
and tail guard of 5.58e-10. The prior-tail contribution is bounded at about
6.03e-22. Low-order quadrature was not accepted because it executed successfully.
Input construction, extended-precision recursion and output summation have
separate rounding allowances. These are checked numerical estimates, not formal
libm enclosures or statistical coverage guarantees.

The actual full-block check validates posterior probabilities. Actual-block
MAP enumeration was not run; path correctness is supported by the exhaustive
fixture and all 27 independent point-limit references. The key failed case's
maximum nonzero posterior is far below .95, so the negative core conclusion
does not depend on an unresolved marginal path choice.

## Conditional findings and disposition

The table uses each first-half territory's maximum nonzero-state probability,
not a replacement for the every-bin core requirement. Boundary-crossing bins
remain included under the earlier start-before-10Mb definition.

| Production case | Executed maximum support | Shared-effect maximum | Shared-effect first-half core |
|---|---:|---:|---|
|20261108|0.997886|0.999999857|Present; bins209–216, minimum support0.994100|
|20261109|0.008886|0.885163|Absent|
|20261110|0.000878|0.715700|Absent|

All alternatives and both paths remain in
[the full 27-comparison results](../dev/regional_amplitude_20261109/panel_v1/results.json).
For production109, whole-chromosome and second-half alternatives also lack a
core, with maximum first-half support about .104 and .107. Production110's
corresponding maxima are about .290 and .028. They are not silently replaced
with the most favourable window. The first case already had initial support;
its distinction is still lost later in the unchanged production workflow.

[Inference] Fixed amplitudes suppress support for these measured weak pool
contrasts, but relaxing that restriction in this one coherent regional model
does not supply sufficient evidence under the unchanged rules. This does not
prove that no RNA information is present or that every continuous-effect model
fails. It does not justify a width search on these outcomes, lower support
thresholds or an assertion of biological clone validity.

No whole-workflow candidate was run, so no new clone-recovery result or error
rate is claimed. Production remains 2/5 on the original cases. All three
recovery failures, the failed state-null false-split target and uncertainty
misspecification remain unresolved. No power or 5% specimen guarantee follows
from these 27 dependent comparisons. The next correction needs additional
identifying information or a justified observation/selection mechanism; this
candidate alone does not pass its immediate improvement gate.

## Complete attempt and resource journal

| Job | Attempt | Result | Process elapsed / peak RSS |
|---|---|---|---|
|23128053|Numerical preflight|Pass,0:0|14.80s /119252KiB|
|23128054|27 fixed-pool comparisons|Pass,0:0|15.75s /118972KiB|
|23128055|Initial independent audit|Fail,1:0: actual block exceeds prespecified12-position bound|1.52s /118328KiB|
|23128524|Whole-block quadrature audit|Fail,1:0 at report serialization of longdouble after numerical checks|8.28s /213656KiB|
|23128638|Serialization-corrected whole-block audit|Pass,0:0|5.21s /211312KiB|

The initial dimensional assumption was wrong: the saved count-half comparison
retains all 19 chr12 positions. The audit stopped as instructed. The
[explicit amendment](../dev/regional_amplitude_20261109/AUDIT_AMENDMENT.md)
retains the same whole block and prespecifies quadrature orders and error checks
before the replacement audit; no easier subset/case or relaxed scientific
criterion is substituted. The first failed source/output remain intact.
The second failure is report encoding, repaired by storing binary64 display
values alongside the original extended-precision guard strings. Numerical
comparisons still use the extended value. No inference rerun was needed.

[Job commands](../dev/regional_amplitude_20261109/jobs.json),
[accounting](../dev/regional_amplitude_20261109/accounting_final.tsv),
[process resources](../dev/regional_amplitude_20261109/resource_summary_final.json),
[preflight](../dev/regional_amplitude_20261109/preflight_v1/report.json),
[accepted audit](../dev/regional_amplitude_20261109/audit_v3/report.json) and
[source preservation](../dev/regional_amplitude_20261109/source_integrity_final.json)
are retained. Every job requests cpu,v6&amd,1CPU,8GB,2h; preflight/panel use
innovation and audits use nslab. No arrays or interactive QOS. All three
754-file frozen source inventories are unchanged, with identical model source
across versions. All 75 production/native/test source files remain unchanged.
The validated403-Python/4-Rust artifact remains the same; its full suite was not
rerun for research-only code. Frozen scientific outputs, defaults, dependencies
and shared environments are untouched.
