# Fitted RNA membership confidence does not repair genomic membership

[Observed pattern, not a guarantee] On all 1200 saved evaluation fits, high
plug-in RNA coassignment confidence remains compatible with substantial genomic
error. For the full-membership arm, mean per-specimen error among pairs assigned
at least95% confidence is38.94% under the one-population state null and42.60%
under the two-clone condition. A95% confidence cutoff therefore lacks support
as the proposed membership correction. This does not establish that every
possible abstention procedure fails.

## Estimand and assumptions

The target is whether two distinct uniformly sampled cells share their generating
CNA population. For frozen RNA responsibilities R, compute the plug-in score
q(i,j)=sum_k R(i,k)R(j,k), including every component and reference-like cell.
This is a conditional RNA-component coassignment score. Interpreting it as genomic
membership requires identified genomic components, conditional label independence
and adequate treatment of fitted-parameter and adaptive-selection uncertainty.
Those assumptions are not established by a large score.

The strongest alternative is confidently separated expression state within one
CNA population. The generating labels distinguish these explanations in this
simulation. Independent corresponding DNA would distinguish them biologically.
No same-cell DNA is added by this experiment; biological n remains zero.

Both frozen joint arms and all three conditions are retained:200 matched-NB null,
200 state null and200 two-clone specimens per arm. These are previously exposed
conditional simulations, not an independent evaluation. The replication unit is
a simulated specimen. Its179700 unordered cell pairs are dependent observations,
not179700 biological replicates.

## Results before the specimen-level selection gate

The diagnostic declares same-population at q>=.95 and different-population at
q<=.05; all other pairs abstain. These are fixed descriptive levels, not a newly
validated clone decision rule or changed scientific criterion. All cells stay in
the denominator. Conditional errors average per specimen; unconditional errors
and coverage preserve abstention. No confidence intervals or significance claims
are made from these exposed results.

| Arm / condition | Pair coverage | Conditional error among retained pairs | Unconditional error | Brier score |
|---|---:|---:|---:|---:|
| Proposal / NB null |93.35%|0%|0%|0.002707|
| Full / NB null |93.35%|0%|0%|0.002707|
| Proposal / state null |42.38%|38.29%|17.09%|0.324559|
| Full / state null |42.06%|38.94%|16.98%|0.325782|
| Proposal / two clones |49.11%|42.59%|20.91%|0.372130|
| Full / two clones |49.10%|42.60%|20.92%|0.372094|

All200 state-null specimens in each arm contain at least one confidently separated
pair. **This is not a specimen false-split rate**: these are unfiltered pairs
before the existing calibrated specimen gate, eligible-group and support rules.
The original specimen error/power event counts and confidence bounds are unchanged.

Brier score is the mean squared probability error. The fixed uninformative score
q=.5 has Brier .25 in every condition, better than either arm under the state-null
and two-clone laws. The fixed one-population score q=1 has Brier0 under either
null but0.500835 under the two-clone law. That baseline illustrates why suppressing
all distinctions cannot qualify as successful recovery. Neither baseline uses
true labels to choose a prediction. These pair metrics are not substitutes for
profile accuracy, cell-level matching, clone detection or the5% specimen criterion.

## Independent numerical and provenance checks

The design was frozen before execution in
`dev/clone_membership_confidence_20261109/DESIGN.md`. The panel verifies chunk
journals and all1200 saved assignment hashes, seed keys, complete600-cell rows,
truth correspondence, and MAP/responsibility consistency. All responsibilities
are finite and nonnegative with row sums within a storage-precision accumulation
bound. Explicit binary64 row normalization changes values by at most4.44e-16;
it is a probability representation conversion, not a model refit. Scores are
never clipped and cells are never dropped.

The independent audit uses matrix multiplication instead of pairwise einsum,
recounts all1200 cases and12000 reliability strata, and verifies summary counts.
Six pair scores agree exactly with rational arithmetic. All144 cell/component
permutation fixtures pass; five invalid-probability fixtures are rejected.
Twelve full-pair scalar summations independently verify Brier arithmetic.
Maximum Brier discrepancy is1.11e-16, within an operation-count binary64 bound.

All execution used compute nodes and the validated noneditable wheel from
`dev/bulk_posterior_export_20261109/environment`; inference was not rerun.

| Job | Account/QOS | CPUs / requested memory | Elapsed | Process peak RSS | Exit |
|---|---|---|---:|---:|---:|
|23071393|innovation|2 /4GB|37s|385348KiB|0|
|23071548|nslab|2 /4GB|13s|95568KiB|0|

Panel evidence: `panel_v1/{cases.tsv,reliability.tsv,summary.tsv,report.json,
assignments.json,run_manifest.json}`. Independent evidence: `audit_v1/report.json`.
Execution inventories and resource records remain in `analysis_JOB/`. Frozen
source, predictions, former evaluation data and scientific criteria are preserved.
Package/native/README/tests remain unchanged from this turn's source snapshot;
only diagnostic code and scientific documentation were added/updated. No full
package suite was rerun because production code was unchanged.

## Decision and remaining work

Do not promote a confidence-only abstention or merging correction. Conditioning
more strongly on fitted RNA labels can preserve confident state subdivisions
without identifying genomic groups. The test rejects that proposed mechanism at
the declared95% operating level; it does not prove nonidentifiability under every
possible independent measurement or observation model. A score-rescaling fit to
these exposed outcomes would not independently validate a corrected procedure.

The original three recovery failures, failed false-split/power tradeoff and
interval-coverage failures remain. Fixing confidence requires independently
supported constraints on nuisance variation or corresponding DNA, plus complete
workflow calibration including selection and fitting uncertainty. Biological
replication and truth resolution remain requirements, not quantities supplied
by additional pair counts or simulations.
