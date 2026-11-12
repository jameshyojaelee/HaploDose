# Gene-level event proposals do not resolve the remaining recovery failure

Two implemented proposal replacements, an NB bin variance derived from gene
moments and a full gene NB likelihood, both recover4/5 original cases. Neither
recovers20261109. Reject both as sufficient recovery corrections and do not launch
their planned larger null/power continuation. Production remains2/5; the earlier
4/5 research procedure remains rejected by its13/100 state-null errors. No default,
test, posterior threshold, false-split criterion or frozen result changes.

## Scientific question and exact intervention

The estimand is distinct regional RNA-profile populations under conditional NB
working models, not genomic clone identity. The question was whether information
lost by bin summation, or the corresponding fitted variance, explains the remaining
failure. Independent gene means within a bin are allocated from native expected
bin means by fitted within-bin gene fractions. This assumes approximately common
gene composition within the parent after the native bin-level nuisance correction.
Depth-dependent composition, reference mismatch, transcriptional state and fitted
nuisance remain competing explanations. Corresponding independently processed DNA
and biological replication would establish the genomic interpretation; n=0 here.

`gene_event_proposal.py` estimates each parent's within-bin fractions from its
actual proposal gene counts only. Reference-only depth-quintile gene dispersion
and all stochastic gene splits remain fixed. The moment arm uses
R_bin=1/sum(fraction_gene²/r_gene), with r_gene already scaled by proposal share.
The gene arm passes those actual gene counts, allocated expectations and sizes
to the existing converged native two-component NB optimizer. Both use the same
native event regions, rankings, candidate deduplication limits and initialization.

All-zero genes have fitted rate zero. A completely unobserved bin makes altered
gene composition unresolved; proposals abstain without a positive floor or oracle
fallback. No source truth, true dispersion, event label or DNA enters fitting.
Reference noise is estimated, so independent count halves are not guaranteed.
Full-count nuisance/parent fitting is retained; proposal-only fraction estimation
does not make the complete workflow independent confirmation.

Weighted confirmation/merging, full-count profile fitting, effective-count gene
veto, posterior.95, credible-core, minimum support and stopping all remain unchanged.
This is a model intervention, not a demonstrated historical arithmetic defect.

## Numerical and exact control evidence

Job23139977 passes exact Fraction moments, count/expectation/variance conservation,
gene permutation, unsupported-bin abstention and direct independent SciPy NB
mixture calculations. Maximum likelihood discrepancy5.83e-16 lies within the
operation/magnitude budget5.97e-12. Two gene allocations with identical count sums
have independently calculated geometric class-likelihood ratios30/49 and24/35;
the native gene responsibilities reproduce them. This proves that summation can
lose information, not that an estimated gene model improves clone discovery.

Control23139978 exactly reproduces every prior numerical array/table in all five
fits and every array in all saved gene-fold files, then passes the unchanged4/5
verdict. This validates the research wrapper before altered arms execute. The sole
existing helper change yields its registry to an explicit research context; old
callers do not consume it. No production package/native/test code changes.

All executed region arrays, selected gene/bin maps, fitted parent compositions,
actual responsibilities, parameters and source matrix paths are retained under
each arm's events/ directory and event_trace.json, linked to exact gene folds.
The four records input → executed source/environment → output → truth/scorer are
available through the new and nested manifests. [Frozen design](../dev/gene_event_proposals_20261109/DESIGN.md).

## Complete original recovery results

| Seed | Original bin control | Moment-bin likelihood | Gene likelihood | Gene assignment accuracy |
|---|---|---|---|---:|
|20261031|pass|pass|pass|.836032|
|20261107|pass|pass|pass|.738095|
|20261108|pass|pass|pass|.756917|
|20261109|fail|fail|fail|.602000|
|20261110|pass|pass|pass|.811623|

All30 criteria remain in each report, including separate normal FPR, tumour TPR,
one S2 group, exactly two S1 groups, original assignment>.65 and chr12-only
differences. The gene arm's110 detects499/510 tumour cells versus500/510 control,
with assignment.811623 versus.814000. It does not improve the existing success.
S2 and S3 each retain one tumour group in every arm.

For109's leading first-half chr12 proposal, mean membership margin is.658452
control,.641899 moment-bin and.672542 gene. All three have zero initial weighted
HMM segments. The gene likelihood does change fitted proposal evidence; its
slightly stronger margin is insufficient to clear the unchanged confirmation
criterion. A merge-only alteration cannot rescue this executed comparison.
There is no justification to lower posterior support or declare unresolved pools
equivalent to genomic clones.

Moment-bin and gene arms fit1701 and1702 event mixtures, respectively. The gain
in representation does not produce the required5/5 recovery improvement under
these estimated nuisance models. This rejects the proposed sufficient correction;
it does not prove that no gene model can help, or that genomic biology is
identifiable from RNA alone. Their fresh false-split/power rates were not estimated
because the prespecified recovery prerequisite failed. The5% simultaneous-bound
requirement is preserved, not passed by skipping those experiments.

## Execution and preservation

Jobs23139977–23139980 all complete0:0, CPU on innovation/nslab, no arrays or
interactive QOS. Numerical gate requests1CPU/4GB/1h; control and altered arms
4CPU/16GB/3h. Scheduler runtimes25s,1m33s,1m21s and3m30s, respectively; process
resources and sampled scheduler memory are both retained. Short-job zero RSS is
missing sampled peak information, not zero memory consumption.

The complete dirty baseline was frozen before editing. Immutable source_v1 runs
against noneditable wheel23138003 read-only (412 Python/4 Rust passes,1 GPU-only
CPU skip). Full package/GPU suites were not rerun for research-only edits. Native
source, existing tests, frozen results and prior decisions remain unchanged.
[All commands, parameters, resource records and outputs](../dev/gene_event_proposals_20261109).
