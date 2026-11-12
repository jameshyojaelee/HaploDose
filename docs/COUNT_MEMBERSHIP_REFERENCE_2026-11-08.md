# Count noise explains little of the saved membership gap

The known-generator count reference reaches67.26% all-cell correctness in the
200 saved two-population simulations, versus68.70% with noiseless regional
amplitudes. Current full-merge assignments reach46.21%; even truth-informed
coarsening of their fixed groups reaches64.24%. These descriptive measurements
separate count noise from partition losses. They do not establish a deployable
correction, biological clone identity or a passed false-split/power gate.

## Scientific contract

The estimand is mean oriented membership correctness over600 cells in each of
200 simulated specimens from the unchanged NB/state generator. This is not
biological replication. Its two populations differ in regional RNA rate by
0.6log(1.5), with Normal(0,.25^2) cell state in the same region and whole-genome
rate normalization. Counts in each of four genomic quarters have equal NB
probabilities conditional on state and depth. Their sums are sufficient, with
NB shape60/.15. The reference integrates the state and uses equal marginal
class priors. It additionally knows the true cell depths, event coordinates,
component effects, reference, dispersion and nuisance distribution.

These unavailable inputs are supplied only to the diagnostic. The comparison
with HaploDose does not isolate a numerical optimizer defect or establish that
the whole gap is recoverable from RNA without independent constraints. The
strongest competing biological explanation remains transcriptional state: a
bimodal programme can generate the same RNA distribution. Independent DNA or
independently justified nuisance restrictions would distinguish genomic identity;
a larger likelihood or more accurate RNA component assignment cannot do so.

## Paired descriptive results

| Quantity | Mean correctness over all cells |
|---|---:|
| Noiseless amplitude, known mixture parameters |68.7033%|
| Observed counts, known true depths/model |67.2642%|
| Same count reference plus known balanced class sizes |67.3267%|
| Joint proposal-merge assignments, one-to-one scoring |46.3917%|
| Joint full-merge assignments, one-to-one scoring |46.2117%|
| Proposal-merge fixed groups, truth-informed coarsening |64.2333%|
| Full-merge fixed groups, truth-informed coarsening |64.2383%|

All cells remain in denominators, including reference-like and unassigned cells.
Extra groups receive no extra credit in one-to-one scoring. Truth-informed
coarsening may map multiple groups to a known true class; it is an attribution
tool, not an executable method. Balanced sorting additionally knows300 cells per
class. The latent and count oracles are empirical references, not per-realization
bounds. These are exposed simulation summaries without new hypothesis tests,
confidence intervals, promotion thresholds or biological validation claims.

The count/noiseless difference is1.44percentage points. The one-to-one/coarsening
gap is18.03points for the full-merge arm. Even coarsening leaves3.03points below
the count reference, under its additional assumptions. [Inference] Further depth
or numerical confirmation alone is poorly motivated for this particular failure;
partition formation and the information available to constrain it matter more.
The previous observational-overlap and specimen-testing limits remain in force.

## Independent computational evidence

Every regenerated count, reference, dispersion, truth, reset, state and mean
array matches both arms' recorded hashes in all200 cases. Quarter aggregation
is checked against full240-bin SciPy NB PMFs across two cells and six state/class
configurations each: the class-independent remainder varies by at most2.28e-13,
within an accumulation allowance derived from operation count and PMF magnitude.

Gauss-Hermite orders16/32/64/128/256 retain the same final membership decisions
at128 and256 across120000 cells. Maximum log-odds changes over successive orders
are0.202854,0.0453503,0.000274799 and2.15e-8. Thus low-order quadrature is not
assumed adequate merely because it executes. Nine selected cell integrals from
cases0/99/199 are independently computed with adaptive integration of full NB
PMFs. Maximum log-odds discrepancy is7.50e-13. Tail bounds use likelihood<=1
and the explicitly omitted standardized Gaussian tails beyond+-12. These checks
support numerical convergence, not a formal global error enclosure. The stated
integration tolerances are numerical requests, not scientific cutoffs.

Independent audit23002207 recounts all200 cases and400 assignment matchings by
exhaustive enumeration. It verifies posterior/balanced decisions, noiseless
reference, source/output hashes, summary means and reported numerical diagnostics.

## Execution, decisions and next constraint

Job23002122 failed exact regeneration onAMDv6 before scoring. No hash criterion
was relaxed. The matchedIntelv4 run23002174 reproduces all200 cases and completes
in1m02s, with947008KiB process peak RSS. The audit completes in3s with83104KiB.
Matching hardware restores regeneration; this is not an isolated causal test of
the architecture difference. Failed and successful runs remain preserved.

No production code, original test, scientific threshold, frozen fit, specimen
gate or comparator changes. The original thinning result remains2/5; its three
failures are unresolved. The broader joint method remains unpromoted after its
failed error/power experiment. Biological nulln=0. Stronger calibration cannot
manufacture membership information or distinguish equivalent biological laws.

Before another membership correction, specify which independently available
observation constrains the currently oracle-supplied quantities. A proposed RNA
fit must remove those oracle inputs and repeat the full observed/null workflow,
including state misspecification and simultaneous power, before any promotion.
Additional allele evidence requires an appropriate ASE/mapping/phase observation
model and matching positive-evidence/null searches; it cannot be presumed genomic.

Implementation is `benchmarks/hardening/count_membership_reference.py` and
`count_membership_audit.py`; locked design, source snapshots, hashes, original
failure, successful outputs and execution/resource records are under
`dev/count_membership_reference_20261108/{reference_v1,reference_v2,audit_v1}`.
The validated wheel is unchanged; the full package suite was not rerun.
