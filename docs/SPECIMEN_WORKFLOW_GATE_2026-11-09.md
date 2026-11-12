# A complete-workflow empirical gate loses the recovery contrasts

The empirical specimen gate is implemented and its calibration banks are built
from344 complete native workflows. It rejects all four previously recovered
original contrasts, leaving0/5 recovery. Reject this candidate as a sufficient
performance correction; the prespecified larger validation is not run. All30
raw recovery values reproduce exactly before gating. No tests, thresholds or
frozen results are weakened to rescue the gate.

## Fixed statistic and calibration argument

The estimand is support for an additional regional RNA-profile distinction within
a specimen, conditional on an explicitly specified null population. It is not
DNA clone identity. Use T=max over final tumour-profile pairs of sum(abs(z))
over that pair's reported regional contrasts. Exclude normal-like profiles by
the native field. No extra profile or no contrast givesT=0. This is a ranking
statistic, not a nominal normal-test p-value. It uses no true labels, event region,
prevalence or DNA. Every null dataset repeats the full public pipeline, including
reference/nuisance estimation, proposals, confirmation, merging and stopping.

Two procedures run on identical RNA: production and the frozen reference-only
depth-quintile weighted/effective-gene candidate. Fresh phase4 has86 replicate
base seeds, matched and continuous-state null conditions. There are172 condition
datasets/344 paired workflows. Replicates are independent within a condition;
conditions and procedures are paired, not independent biological specimens.

Each procedure's additional distinction must exceed the largestT in both null
banks. Equivalently, both inclusive rank p-values must be at most1/87; ties do
not pass. Derive the bank size exactly: for n iid null scores, the probability
that the maximum leaves more than5% upper-tail mass is at most.95^n. There are
two procedures × two nulls, so require.95^n<=.05/4. The smallest n is86;85 fails.
The exact Fraction checks pass, including tied, all-zero and strict-exceedance
fixtures. A union bound gives failure probability at most.0485613 across those
four calibration tolerances; independence between banks is unnecessary.

This statement assumes the observed score's null law matches the calibrated law.
The native synthetic RNA generator, nuisance distribution and sampling regime
are explicit assumptions. It is not biological error control, empirical held-out
validation or control of unsupported subdivisions inside an already multiclonal
tumour. Allele/spatial configurations are not tested. Continuous states from
another law, reference mismatch or different composition may invalidate transport.
Independent corresponding DNA and biological null replication remain necessary.

[Fixed design](../dev/specimen_workflow_gate_20261109/DESIGN.md),
[statistic, bank construction and exact checks](../benchmarks/hardening/specimen_workflow_gate.py)
and [recovery application](../benchmarks/hardening/specimen_workflow_gate_recovery.py)
retain the implemented procedure.

## Thresholds and complete original recovery probe

Matched-null maxima are0 for both procedures. State-null maxima are18.285016
production and20.906054 candidate. The inference runner exposes the same statistic
in observed and simulated outputs, with all source/helper/environment identities
recorded. Source truth metadata are excluded from inference inputs.

The five original inputs are rerun on current immutable source. Every raw
criterion and original assertion verdict matches the earlier4/5 candidate.
Apply the frozen candidate bank to the new outputs, preserving tumour calls
but marking additional structure unresolved when the gate fails. The separate
view uses one unresolved RNA population, not a certified single genomic clone.
Raw model outputs remain intact; original scoring and assertions are unchanged.

| Seed |NativeT|State-bank inclusive rank p|Gated original recovery|
|---|---:|---:|---|
|20261031|18.948355|.068966|fail|
|20261107|18.036617|.080460|fail|
|20261108|16.409896|.103448|fail|
|20261109|0|1|fail|
|20261110|16.070396|.114943|fail|

All four raw successes lose their required two-group/assignment/chr12 criteria.
The109 case already lacked a qualifying raw distinction. Even a nominal.05
rank gate on these same banks would reject the four contrasts; this observation
does not authorize changing the locked threshold. No independent validation
false-split or power rate is estimated for the rejected gate. Its planned phase5
pilot/final400 expansion is **not run** because the recovery prerequisite fails.

An important transport limit remains: these calibration queries contain600
tumour cells with no normal admixture, whereas the original recovery queries
contain510 tumour/90 normal cells and different specimen configuration. The
recovery probe therefore tests transport as well as the gate. A conditional
calibration statement for its own null laws does not prove calibration for the
original mixed specimens. This is a justified next nuisance-matching question,
not proof that all possible RNA statistics lack power or that biology is
unidentifiable in every assay. Do not retune a threshold using these outcomes.

[Bank report](../dev/specimen_workflow_gate_20261109/bank_v1/report.json),
[all30 raw/gated recovery criteria](../dev/specimen_workflow_gate_20261109/recovery_v1/report.json)
and each separate unresolved view preserve the evidence. These simulation probes
provide no new biological replicate or DNA truth.

## Replication-label correction

The preceding phase3 summary's frozen `independent_simulations=300`(pilot30)
overstates independence across conditions. The generator uses the same base seed
and reference counts for matched/state/two-clone conditions within each replicate:
for example all three phase3 replicate0 records have seed202611093000. Correct
counts are300 condition datasets/100 independent base seeds, and pilot30/10.
Each condition still has100 independent replicates; its exact binomial bounds
and Bonferroni decisions are unchanged. Bonferroni does not require independent
conditions. Original JSON/results remain frozen. Future summary source fields
now say `condition_datasets` and `independent_base_seeds`; the earlier report has
an explicit amendment rather than silently replacing historical evidence.

## Execution and decision

Eleven CPU jobs complete0:0: calibration23134196–23134203, bank23134204,
raw recovery23134205, gated recovery23134206. Eight paired calibration chunks
use4CPU/16GB/2h, four innovation/four nslab; raw recovery innovation4CPU/16GB/2h;
bank/application nslab1CPU/4GB/1h with explicit dependencies. Process peak
2120620KiB. Immutable source_v1, validated noneditable wheel and neutral-directory
imports; no arrays, interactive QOS, new dependencies or frozen-output writes.
[Commands/dependencies](../dev/specimen_workflow_gate_20261109/jobs.json),
[accounting](../dev/specimen_workflow_gate_20261109/accounting.tsv),
[actual resources](../dev/specimen_workflow_gate_20261109/resource_summary.json).

All75 production/native/test files remain unchanged. Full package/GPU suites
were not rerun for research-only code. No package default or scientific criterion
changes. Keep the failed candidate, exact calibration argument, transport limit
and all recovery failures visible. The original three production failures and
biological5% false-split/power objective remain unresolved.
