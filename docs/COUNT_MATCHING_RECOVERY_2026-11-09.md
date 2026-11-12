# Two depth-matching corrections fail to rescue the remaining confirmation case

Exact averaging of the existing count-matching kernel and an alternative physical
read-depth kernel both fail the locked case109 continuation gate. Its maximum
nonzero state support increases from.000954 to.003890 or.004685, respectively;
neither produces a qualifying core. All67 initial comparisons reproduce exactly.
Do not integrate either candidate or spend a larger full-workflow/null/power batch
on them. Production recovery remains2/5; the earlier4/5 research candidate remains
rejected for excess expression-state false splits. No criteria or frozen changes.

## Estimand, observation law and competing explanations

The target is a distinguishable regional relative RNA effect between proposed
cell populations, conditional on the existing discovery model. It is not DNA
copy number or independently validated clone membership. The actual initial
chr12 proposal in seed20261109 has strong working support before depth matching,
but loses it after matching despite retaining all502 positive-weight cell IDs.
The question here is whether ancillary thinning or its observation law explains
that loss. No confirmation outcome selects a random redraw, weight or new effect.

First candidate: for fixed NB(E,r), the existing matching kernel has
K|Y~BetaBinomial(Y,wr,(1-w)r), giving marginal NB(wE,wr). Average the exact
alternative/null likelihood ratio over every integer K=0..Y, using a normalized
log-probability recurrence. Conditional averaging uses every observed cell and
gene while respecting the EXACT executed matching weight w. It does not average
log likelihoods or posteriors, bypass matching or accept initial evidence alone.

Second candidate: physical read/molecule downsampling K|Y~Binomial(Y,w) gives
NB(wE,r), retaining the original cell-level Gamma-Poisson dispersion. This is
a different experiment from an independent Gamma split. The beta-binomial law
is not universally incorrect. Discovery/confirmation fold splitting stays
unchanged; only auxiliary matching changes. The discarded matching counts are
not used in another confirmation test, so retained/discarded dependence is
permitted conditionally. Exact conditional LR averaging has a closed form.

Both retain the same component-normalization proportions, discovery membership
probabilities, selected genes/cells, expected counts, priors, genomic boundaries,
every-bin.95 support and minimum-five-bin rule. Weight1 follows the original gene
calculation with its exact operation order; weight0 contributes no evidence.
Neither adds a noise constant or selects among draws. Cell/gene vetoes, final
merge rules and the5% specimen false-split/power objective are untouched.

For FIXED correctly specified nuisance, weights and proposals independent of
confirmation observations, the averaged LR has null expectation1 by iterated
expectation. Actual full-data depth matching, nuisance estimation, dependence
and proposal selection prevent importing that identity as a workflow guarantee.
Normalized LR state weights are working support, not generative posterior
probabilities for unthinned observations. Regional transcriptional state,
reference mismatch and count-model misspecification remain biological/statistical
alternatives. Corresponding independent DNA at the claimed resolution and
biological replication would establish genomic validity; biological n=0 here.

## Decisive design and independent numerical references

The [averaging design](../dev/averaged_count_matching_20261109/DESIGN.md) and
[physical-depth design](../dev/binomial_depth_matching_20261109/DESIGN.md) precede
their method results. Both require an actual depth-matched109 chr12 qualifying
core before any full recovery run. A pass would still require the separately
announced historical exact control, all30 original recovery criteria and fresh
complete-workflow error/power validation. The failed prerequisite stops both.
Earlier failed historical controls are not replaced or reinterpreted.

Independent Fraction likelihood/kernel fixtures, Decimal60/90-digit convergence,
and complete SciPy enumeration of all retained counts, cell allocations and state
paths verify the beta-binomial implementation. A rational enclosure bounds the
fixed-null omitted expectation by4.94e-22. Maximum numerical discrepancy is
7.11e-15. Input rejection, exact endpoints and cell/gene permutation also pass.

Physical-depth checks independently enclose exact NB marginal probabilities
4/9,8/27,4/27 with an omitted bound3.44e-71, and verify retained/discarded
covariance1/2 in the specified fixture. Complete SciPy enumeration and separate
extended-precision scalar formulas agree, maximum discrepancy3.56e-15; actual
extended mantissa is63 bits. Endpoint and permutation checks pass. Allowances
derive from operation counts, scales, finite-tail bounds and independent
reference discrepancies, not a universal epsilon or a scientific tolerance.

Implementations:
[averaged_count_matching.py](../benchmarks/hardening/averaged_count_matching.py),
[binomial_depth_matching.py](../benchmarks/hardening/binomial_depth_matching.py).

## Exact saved-data replay and every result

Reconstruct original confirmation counts from the15 saved source/prop CSR folds
using literal global cell and gene identities. Verify initial observed counts,
actual thinned-count conservation, matching weights, expected-count scaling and
source NB sizes. All180 original trace entries remain linked:88 gene-confirmation,
88 ordinary confirmation and4 merge entries. Preserve21 gene abstention/fallback
entries and the92 other-stage records. The67 identified initial comparisons plus
five actual matched comparisons give72 replayed comparisons in each candidate.
There are185 output records because identified entries can contain both contexts.

| Actual109 matched chr12 comparison | Original draw | Averaged beta-binomial | Averaged physical-depth |
|---|---:|---:|---:|
| Maximum nonzero support |.0009540565|.0038895411|.0046850102|
| Qualifying core |None|None|None|
| Represented cells |502|502|502|

The source event is bins209–217, maximum original gene count19, complete finite
kernel support60239 terms across those observations. All67 initial posterior
arrays are exact in each experiment, maximum discrepancy0. No matched comparison
changes whether it has a core. This falsifies either intervention as a sufficient
recovery correction; it does not establish that every possible depth model is
irrelevant. No new full-workflow recovery, merge, false-split, power or biological
uncertainty result is reported. Those checks are explicitly not run after failure.

[Averaged results](../dev/averaged_count_matching_20261109/study_v2/report.json),
[physical-depth results](../dev/binomial_depth_matching_20261109/study_v1/report.json).

## Compute, failed attempt and preservation

| Job | Experiment | Account/QOS | CPUs / RAM / requested time | Outcome |
|---|---|---|---|---|
|23150637|beta-binomial numerics|nslab|1 /4GB /1h|complete0:0,24s|
|23150638|first saved-case reader|innovation|1 /8GB /1h|failed on heterogeneous trace stage|
|23150650|corrected saved-case reader|innovation|1 /8GB /1h|complete0:0,69s; scientific gate fails|
|23150756|physical-depth numerics|nslab|1 /4GB /1h|complete0:0,10s|
|23150757|physical-depth saved cases|innovation|1 /8GB /1h|complete0:0,67s; scientific gate fails|

The first reader incorrectly expected a case_seed on ordinary confirmation trace
entries. It stops with KeyError before a complete scientific report. Source_v2
recognizes those original stages explicitly; it does not remove any gene case or
relax a criterion. Failed output and both snapshots stay intact. There is no
historical control/node-retry search or random redraw selection.

Every executable check runs on CPU compute nodes using the existing read-only
noneditable validated wheel, native identity, neutral directory and recorded
backend/allocation/resources. Peak process RSS and accounting are saved per
job; scheduler RSS0 is not a substitute for the process measurement. Job commands,
source manifests and immutable snapshots reside in each task directory. Shared
installations, production/native source, original tests and frozen outputs stay
unchanged. The prior421 Python/4 Rust pass and one GPU-only skip are separate
evidence; no full package suite or GPU parity rerun is claimed for research code.
