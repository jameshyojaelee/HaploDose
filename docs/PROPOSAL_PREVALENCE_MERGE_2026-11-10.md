# Test whether imported mixture proportions cause insufficient merge support

COMPLETE saved-output gate: reader23170951 exits0:0 in28s,with process
peak231044KB. ALL4arrays/28tables including logs and1fold archive are EXACT;
both original effective-support records,all scalar prevalence fits and both
saved responsibility matrices reproduce EXACTLY. Input hashes remain unchanged.

Reject proposal-fold prevalence profiling as a sufficient correction. On the
two existing merge attempts,the smaller conditional tumour prior falls from
.096419/.096423 to.076747/.076749. MAP minority sizes fall35→30 and effective
supports fall25.065572→21.893102 and25.066589→21.893753. Neither pair is
testable at the unchanged30minimum,and neither has a pre-veto segment. This
falsifies the proposed sufficient-correction mechanism on this case. Do not
implement it as a merge rule,retune a threshold or expand to a prevalence sweep.
Case2 remains unresolved; this negative analysis adds no recovery or error claim.

Execution failure retained:23169186 exitsFAILED1:0 after23m34s. Both nested
inference arms and their pilot manifest are COMPLETE,but the final checker
tries to import refinement_checked_replay.py,omitted from source_v1. No method
or scientific check failed there;the required output gate had not executed.
The new standalone reader runs the FULL unchanged checks from saved outputs,
without rerunning inference or overwriting the failed parent manifest. Its
source_v2 preserves all913executed files exactly and adds only reader/design.
[Complete gate and negative result](../dev/proposal_prevalence_merge_20261110/summary_v1/report.json),
[reader manifest](../dev/proposal_prevalence_merge_20261110/summary_v1/run_manifest.json),
[failed checker preserved](../dev/proposal_prevalence_merge_20261110/diagnostic_v1/run_manifest.json),
[reader design](../dev/proposal_prevalence_merge_20261110/READER_DESIGN.md).
Earlier preparation checkpoints below retain their original context.

One observation-only conditional diagnostic23169186 is RUNNING, CPU/innovation,
4CPU32GB90h on original005/runtime003. It follows the complete exact case2
trace, whose proposal-fold merge assignments565/35 leave25.066 effective
minority cells, below the unchanged minimum30. No inference output is changed.

Condition on the actual proposal-fold likelihoods and fitted component rates/
nuisance. Retain the fitted normal mass and maximize only the relative mixing
proportion of the two tumour components. This is a concave scalar likelihood;
the existing adjacent-float profile includes boundaries and abstains when flat.
Compare original versus profiled responsibilities under the SAME MAP eligibility,
eight depth strata, classification margins, squared-weight sampling terms and
30-cell requirement. Save every likelihood, row ID, prior, responsibility and
pre-veto segment. No added gene vetoes or scientific decisions are committed.

Numerical checks COMPLETE on the compute node: the independent rational1/8
maximizer and SciPy derivative root agree with zero discrepancy. Component
exchange, exact boundaries, flat abstention and additive cell-wise likelihood
constants pass. The operation-count arithmetic allowance2.55795e-13 applies to
these identities, not to a scientific support threshold.

Interpretation awaits ALL4arrays/28tables/logs/folds matching the archived
workflow EXACTLY and original effective support being reproduced exactly. If
profiling leaves support insufficient or removes the pre-veto contrast, reject
this explanation as a sufficient correction. A likelihood increase alone is
guaranteed by fitting and does not establish benefit. Passing pre-veto support
would still require unchanged gene vetoes and complete recovery/error/power
testing before adoption. Full-count nuisance remains fitted; no independence,
genomic-identification, calibration or biological claim. Biological n=0.

All911prior immutable source files remain exact in the913-file diagnostic
snapshot. The old validated environment is read only; no package/default/frozen
result changes. [Exact design](../dev/proposal_prevalence_merge_20261110/DESIGN.md),
[implementation](../benchmarks/hardening/proposal_prevalence_merge.py),
[numerical result](../dev/proposal_prevalence_merge_20261110/diagnostic_v1/numerics.json),
[job command](../dev/proposal_prevalence_merge_20261110/jobs.json).
