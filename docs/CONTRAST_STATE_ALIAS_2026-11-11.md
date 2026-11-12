# Why the contrast-only model loses state-free profiles

The complete contrast-only ablation fails28/32 RNA and52/64 simultaneous tests.
In clean high-depth data,RNA MSE improves to.000270567 while state-free MSE
worsens to.000919615,versus original.000340596/.000167622. This diagnostic tests
whether the fixed pair constraint moves missing contrast into continuous state.
It cannot rescue the rejected method or supply a deployable oracle.

Use all48 original clean-case seeds at both depths,two offsets for contrast-only
and coupled models,and the original unrestricted model:480 saved fits. Reconstruct
responsibilities and state means independently from full multinomial probabilities.
Verify the algebraic decomposition of expected centered log probabilities into
shared mean,discrete contrast and continuous expression state. Compare the known
contrast with the pair-constant space and inspect both offset errors. Exact labels
are diagnostic-only;all outputs were locked before their use. Biological n=0.

If the incompatible offset loses discrete contrast while its continuous term
carries the corresponding signal,that supports allocation as the cause. Similar
aligned-offset failure or absent state contribution rejects that as sufficient.
No formal new significance family,refit,penalty,threshold or selected comparator.

Job23220208 requests one nslab CPU,4GB,one hour. Public code,all original criteria,
frozen results and342 prior ledger values remain. Genomic claims still require
independent DNA;model decomposition does not identify biological state.

[Frozen design](../dev/contrast_state_alias_20261111/DESIGN.md),
[executed trace](../dev/contrast_state_alias_20261111/trace.py),
[source identities](../dev/contrast_state_alias_20261111/SNAPSHOT.json),
[submission](../dev/contrast_state_alias_20261111/submission.json),
[rejected full ablation](CONTRAST_ONLY_JOINT_2026-11-11.md).

## Completed mechanism

Job23220208 completes0:0 in45.92s,92032KiB process RSS and73164KiB scheduler
RSS. All480 cached fits verify. The maximum independent state-mean discrepancy
is2.40e-12 and centered log-probability identity discrepancy9.77e-15. No model
was refitted;all342 prior ledger values and91 protected files remain exact.

At depth3000,the contrast-only aligned offset selectsK2 in48/48,with mean
discrete true-contrast projection1.005 and continuous projection0.000503.
The shifted offset selectsK2 in0/48;the continuous projection is0.8894 while
the discrete term iszero. The shifted pair space cannot represent50.54% of
the centered true contrast's squared norm;its continuous term recovers88.30%
of that residual. State-free whole MSE is.000129433 aligned versus.00341834
shifted. Original unrestricted joint selectsK2 in48/48,discrete projection1.007.

This directly supports boundary-induced allocation into state in the saved
working model. It does not prove DNA identity in biological samples. At low
depth,even aligned contrast-only selectsK2 inonly9/48 and the original in0/48,
so correcting the grid alone cannot be assumed to solve all recovery limits.
The general correction must respect data-supported structural uncertainty,
without choosing a grid using truth or interpreting RNA prediction as DNA.

[All480 decompositions and ten summaries](../dev/contrast_state_alias_20261111/run_23220208/results/REPORT.json),
[actual resources](../dev/contrast_state_alias_20261111/run_23220208/resources.txt).
