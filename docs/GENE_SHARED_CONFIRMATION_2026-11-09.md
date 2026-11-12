# Shared membership retains an initial contrast but does not repair recovery

Four implemented exact shared-membership confirmation candidates each recover4/5,
failing20261109. Two integrate bin/gene evidence on an integer depth-balanced cell
subset. Two retain cells and match depth by NB count thinning. All fail the frozen
5/5 continuation prerequisite. Reject them as sufficient corrections; their
larger null/power continuation is not run. Production remains2/5, the earlier4/5
procedure remains rejected for13/100 state-null errors, and no default changes.

## Working model and identifying limits

The estimand is distinguishable regional RNA populations, conditional on fitted
proposal probabilities and nuisance parameters. Each cell has a single latent
binary membership shared across all genes/bins of a candidate event. Conditional
on membership and the same three-state genomic-region path, gene counts have
independent NB working likelihoods. Gene log likelihoods sum within bins before
state enumeration, so gene multiplicity creates no additional transitions.

Within-bin composition is estimated from the actual parent proposal counts only;
gene dispersion is estimated from fit-reference depth quintiles only. Bin sizes
match summed gene NB moments. Exact latent integration does not validate these
nuisance estimates, create independent confirmation, or distinguish genomic
dosage from a continuous regional transcriptional state. Corresponding independent
DNA and suitable biological replication would distinguish that biological claim;
biological specimens added here=0. These original simulations are development
regression cases, not population error estimates.

`gene_shared_confirmation.py` integrates memberships and every state path exactly
up to the fixed one-million-path cap. Longer proposals retain the old procedure
with an explicit fallback record. Each arm has88 event-candidate entries and19
such fallbacks. Moment-bin arms execute74 exact comparisons; gene arms execute72
and abstain twice because a gene unobserved in discovery has confirmation counts.
No positive expression floor fills that missing composition information.
Discontinuous/uncallable regions cannot acquire artificial supporting evidence.

Native proposal regions, rankings, deduplication, full-count initialization and
refinement, weighted merge, effective-cell/gene vetoes, posterior.95, credible-core
length, state effects, prior/transitions and all original recovery assertions stay
fixed. Original pooled contrasts and z scores remain descriptive ranking summaries;
the new model replaces only candidate-confirmation supported segments.

## Decisive recovery and matching experiment

Exact control23143870 reproduces every retained numerical array/table, every gene
fold and all30 original criteria before candidates execute. Every new candidate
has the same original criterion values as that4/5 control. All retain one S2 and
one S3 tumour group in each seed. Seed109 still has one S1 tumour group, original
assignment.602 and no final chr12 difference. Seeds31,107,108,110 retain their
previous passes and assignment.836032,.738095,.756917,.814000, respectively.

The table describes109's actual returned chr12 first-half proposal, bins209–217.
The reported posterior is its largest nonzero-state probability, a diagnostic
bound rather than a replacement for the unchanged every-bin credible-core rule.

| Model | Initial maximum | Initial cores | Matched maximum | Matched cores | Final recovery |
|---|---:|---:|---:|---:|---:|
| Moment-bin, integer cell balance |.999842|1|.082241|0|4/5|
| Gene, integer cell balance |.999899|1|.107874|0|4/5|
| Moment-bin, count balance |.999842|1|.003551|0|4/5|
| Gene, count balance |.999899|1|.000954|0|4/5|

The initial gene core spans209–216 with minimum posterior.999748. Passing this
stage does not establish full recovery. The integer-cell result motivated a
follow-up retaining every positive-weight cell rather than discarding cells.
Both count arms retain all502 proposed cells in this comparison and still fail.
Thus cell removal alone is not a sufficient explanation of lost matched support.
It does not justify bypassing matching, accepting the initial test alone, selecting
a favorable redraw, or lowering support thresholds. Merge-only changes cannot
restore this candidate, which never passes matched confirmation.

The count follow-up uses the SAME executed depth/membership matching weights w.
For NB count Y with size r, draw P~Beta(r*w,r*(1-w)), then Y'~Binomial(Y,P), with
deterministic endpoints. Under the working law the retained count has mean w*E
and size w*r. Each positive-weight cell retains its identity/shared membership;
the original effective-cell veto still applies. Parent component-normalization
proportions remain fixed. SHA256 derives a reproducible seed from the fixed
method tag, case seed, ordered global rows, source and native seed. Actual
thinned counts, weights and seed are saved. There is one prespecified draw;
no redraw search follows this failure.

## Independent calculation and execution

Numerical23143869 explicitly enumerates all9 paths ×8 binary cell allocations
for a three-cell/four-gene/two-bin fixture under all four boundary conditions,
using independent full SciPy NB probabilities. Maximum evidence/posterior/
membership discrepancy is1.08e-15, within the operation/magnitude budget3.31e-11.
Cell/gene permutation and fixed-rate subsetting agree. Executed gene counts and
expectations conserve bin totals, checked with count equality and operation-based
floating-point allowance. Previous exact same-sum gene information fixtures remain
linked; they do not imply a recovery benefit.

Numerical23144894 repeats those checks and adds independent exact Fraction
beta-binomial probabilities for Y=6,r=2,w=1/3. Conditional mean2 and variance32/9
are exact. SciPy probability discrepancy1.67e-16 is within the product-operation
budget8.08e-14. Thinning endpoints, conservation and deterministic repetition pass.
These are checked numerical budgets, not formal bounds on library functions or
statistical confidence intervals.

| Job | Role | CPU/RAM/limit | Scheduler elapsed | Exit |
|---|---|---|---:|---:|
|23143869|joint-model numerics|1/4GB/1h|26s|0:0|
|23143870|exact retained control|4/16GB/3h|1m32s|0:0|
|23143871|moment-bin, integer matching|4/16GB/3h|2m14s|0:0|
|23143872|gene, integer matching|4/16GB/3h|2m09s|0:0|
|23144894|joint-model and count-thinning numerics|1/4GB/1h|22s|0:0|
|23144895|moment-bin, count matching|4/16GB/3h|2m09s|0:0|
|23144896|gene, count matching|4/16GB/3h|2m25s|0:0|

All use CPU on innovation/nslab, without arrays/interactive QOS. Baseline is the
complete dirty pre-edit source. Immutable source_v1 executes the first four jobs;
source_v3 executes the follow-up. Unused source_v2 is retained before explicitly
expressing the probability comparison's operation-based numerical allowance.
There were no failed executable jobs. Source/input/environment/backend/resource
records link actual counts, proposals, selections, likelihoods, posteriors and
memberships to each output and scorer. Missing scheduler RSS on short jobs is
not zero memory use; process peaks are separately recorded.

The installed noneditable wheel23138003 remains unchanged, with prior412 Python/
4 Rust passes and1 GPU-only skip. Full package/GPU suites were not rerun for
research-only files. The5% whole-specimen false-split criterion and simultaneous
power assessment are preserved but remain unmet, not implicitly passed by stopping
these candidates before null experiments. No genomic validation or SOTA claim.

[Design and follow-up frozen before execution](../dev/gene_shared_confirmation_20261109/DESIGN.md),
[all seven commands and dependencies](../dev/gene_shared_confirmation_20261109/jobs.json),
[gene count-balance results](../dev/gene_shared_confirmation_20261109/count_gene_v1/report.json).
