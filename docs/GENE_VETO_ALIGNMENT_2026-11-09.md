# Aligning gene and regional pools restores the last research recovery case

The aligned-pool research candidate passes all five original recovery cases,
versus4/5 for the same known-gene-noise, membership-weighted comparator.
20261109 recovers two tumour groups, assignment accuracy0.779798 and only the
chr12 distinction. Production remains2/5 and is not changed. These exposed
simulation cases establish a bounded recovery mechanism, not genomic accuracy,
specimen error control, calibrated uncertainty or superiority to other packages.

## Primary measurement chain and decisive observation

The fixed task `dev/gene_veto_weights_20261109/DESIGN.md` investigates a
population mismatch before changing inference. The weighted regional comparison
uses depth-matched membership weights, while gene_consistent pools the same
test-half genes with hard labels. A wrapper captures the actual sparse W object
passed from matched_differences to pooled, keyed by specimen, labels and pair.
Its veto returns remain unchanged. Fitted cells/differences and all30 original
criterion values exactly replay the frozen both-weighted arm.

Of81 gene-veto calls,38 have corresponding captured matched pools covering the
gene subset. The only differing decision is20261109's later merge comparison,
bins209–216. Its regional HMM already has a resolved segment. Hard gene pools
give49/79 direction-concordant genes and nominal sign-tail0.021082722, failing
the existing strict0.01 rule. Gene counts pooled with the exact regional weights
give45/60 genes and nominal tail0.0000672570405, passing that rule. Both use
the exact same gene test counts. The other37 matched decisions agree; the
remaining43 calls, mostly shared profile comparisons, have no matching pool and
are explicitly unavailable for this diagnostic.

The alternative explanation is selected/transcriptional separation rather than
a genomic clone. Changing the pool can change the measured contrast without
making the binomial sign tail valid after selection or under gene dependence.
The nominal values above are mechanism diagnostics, not calibrated p-values.

## Fixed intervention and unchanged recovery requirements

`dev/gene_veto_aligned_20261109/DESIGN.md` is frozen after the diagnostic audit
and before the candidate. The research hook substitutes the gene-veto decision
only when the actual corresponding regional pool is represented completely by
the gene evidence. Original early cell-support abstention and unmatched behavior
remain. It retains8-gene support, count>=10 eligibility, .5 pseudocount and
strict alpha0.01. Weighting changes eligible gene counts and the estimand;
inheriting those constants does not validate them on fractional weighted counts.

The model, membership weights, gene-native splits, true source size10, gene test
half, HMM support, minimum bins/cells, profile fitting, nuisance estimation,
search/stopping and original six recovery assertions remain fixed. The former
hard return is still recorded. No candidate is protected from terminal validation.

| Seed | Fixed-pool comparator pass | Aligned candidate pass | Comparator accuracy | Candidate accuracy |
|---|---|---|---:|---:|
|20261031|yes|yes|0.836032|0.836032|
|20261107|yes|yes|0.738095|0.738095|
|20261108|yes|yes|0.790514|0.790514|
|20261109|no|yes|0.601202|0.779798|
|20261110|yes|yes|0.810000|0.810000|

Accuracy is the original conditional-on-tumour-call assignment criterion. In
20261109 normal false calls remain2/90, tumour detection declines499/510→495/510,
and the separate two-group, one-group-S2 and chr12-only criteria all pass.
These are five previously exposed development cases, not independent patients
or prospective performance estimates. No significance claim is made.

The estimand is the direction of regional relative RNA abundance in the same
weighted cell populations as the regional comparison. Raw gene counts obey the
simulator's independent gene NB law before weighting. Interpretation as genomic
dosage or clone identity requires adequate reference/nuisance controls and DNA
evidence at the corresponding resolution. Known dispersion remains an oracle;
fitted full-count nuisance, weights and adaptive selection remain coupled.

## Independent calculations and execution

Diagnostic23119911 and audit23119928 complete; a stronger audit23120056 directly
compares all15 fold arrays with the frozen comparator. Candidate23120037 and
audit23120038 complete, followed by direct-fold audit23120057. All six jobs
exit0 on CPU innovation/nslab, no arrays. Replay/candidate requests4CPU/8GB/2h;
audits1CPU/8GB/1h. Actual maximum process RSS1634220KiB, below1.56GiB. Runtime,
memory, installed wheel/native identity, CPU backend and exact arguments are in
each task's jobs.json, accounting_final.tsv and analysis_JOB records.

Final audits independently reconstruct81 diagnostic and84 candidate veto calls.
Sparse matrix products are checked by columnwise compensated scalar sums, with
nonnegative-sum operation-count precision budgets. Median-centered log-ratio
signs are independently computed by rational odds comparisons: for even-size
medians compare each squared odds with the product of the two middle odds.
Binomial tails use exact integer combination sums, not a floating cutoff.
Original hard decisions and every applied candidate decision agree. Both sets
of30 original recovery criterion values are independently recounted, including
the comparator's scientific failures. Each replay's15 saved gene folds equal
the frozen comparator arrays, including gene veto subsets, mapping and counts.

The diagnostic and candidate retain their first audit outputs. audit_v2 adds
explicit transitive hashes and direct15-fold equality, without rerunning inference.
All76 protected package/native/test/README/build files match the dirty source
baseline. No existing tests or frozen artifacts change. Full package tests are
not rerun for these research-only scripts; the installed wheel's previous
395 Python/4 Rust passes and1GPU-only skip remain separate evidence.

Added research files are benchmarks/hardening/gene_veto_weights.py and
gene_veto_weights_audit.py. Executed source snapshots and output hashes reside
under the two task directories. Observational and applied modes are explicit;
the public package does not acquire a new option or default.

## Decision and next gate

Freeze the candidate and same-fold comparator for a new complete-workflow
null/power experiment. This must reproduce initial grouping, nuisance fitting,
all proposal/confirmation/gene vetoes, merges and stopping. It must include a
matched one-CNA population, a one-CNA population with correlated expression
state, and true CNA contrasts under the same nuisance. Report exact group count,
assignments and abstention alongside simultaneous specimen false-split bounds
and detection power. No thresholds may be chosen from that evaluation's outcomes.

The original5% whole-specimen upper confidence bound and joint power requirement
remain unchanged. No new null/power or repeated-workflow coverage test has run
in these two tasks. Fitted membership/selection uncertainty and estimated gene
dispersion remain unvalidated. Biological n=0: suitable independent null specimens
and corresponding DNA with documented sensitivity remain missing. Earlier
joint-model state-null failures and uncertainty failures are retained. This
candidate is not promoted, and the manuscript/SOTA goal remains unachieved.
