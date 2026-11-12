# Preserve molecule-split information while modelling shared sampling noise

The scientific goal remains accurate RNA-derived genomic profiles and supported
CNA-defined populations, validated independently at the claimed resolution.
This new performance candidate changes the observation split and its likelihood,
instead of changing a confirmation, support or false-split threshold. It uses
binomial molecule allocation and explicitly shared gene sampling noise. The
entire fixed nine-case recovery/state screen and its independent checks are
recorded here. This is research code; no package default is promoted.

## Why this is a distinct intervention

The [linear-depth experiment](LINEAR_DEPTH_STATE_2026-11-11.md) retains four true
contrasts and rejects all four known state controls, but still rejects required
20261109. Its supported minority and numerically resolved score rule out a
support-boundary or rounding explanation for that executed failure. Earlier
gene beta-binomial thinning and reference-only noise candidates remain separate
negative/partial results. This experiment does not repeat their split seeds or
silently change their failed verdicts.

Beta-binomial thinning can yield independent NB parts under the specified full
NB law and matching shape. Plug-in shape, heterogeneous genes and shared state
do not inherit a blanket independence guarantee. Binomial molecule allocation
instead preserves shared gene-level sampling variation. Its pieces are
dependent, so their independent-NB product is the wrong joint model. The new
candidate deliberately models that dependence.

For one gene,Y~NB(mu*exp(eta),r),D|Y~Binomial(Y,f),C=Y−D:

    D ~ NB(f*mu*exp(eta),r)
    P(D,C|eta)=NB(D+C;mu*exp(eta),r)*Binomial(D;D+C,f).

The shape stays r; multiplying it byf would describe a different marginal.
The binomial factor is independent of eta and cancels in the competing-model
ratio, even after integrating the cell state. Conditional evidence is therefore
the full-count NB model difference minus the discovery NB model difference.
It is not a product of separately redrawn latent states or independent gene
halves. Both gene sampling intensity and the shared cell-state dependence are
accounted for under the working law.

The estimand is conditional confirmation-RNA log-predictive separation between
two regional populations and one Gaussian continuous-state population. It is
not DNA dosage, a genomic membership probability or a calibrated specimen test.
The observation model remains gene NB with estimated reference expectations
and noise, plus the frozen linear outside-event depth trend. Gene-specific
noise mismatch, correlated programmes, reference transport and RNA-equivalent
transcriptional states are competing explanations. Independent technical
sampling measurements can test the noise law; same-cell DNA and independent
biological replication can test genomic membership and biological error.
No such new biological observation is added here.

## Frozen design and primary evidence chain

[Design](../dev/shared_gamma_split_20261111/DESIGN.md),
[split implementation](../dev/shared_gamma_split_20261111/split.py),
[shared-noise prediction](../dev/shared_gamma_split_20261111/joint.py),
[unchanged model family](../dev/shared_gamma_split_20261111/frozen_model.py),
[driver](../dev/shared_gamma_split_20261111/study.py),
[source hashes](../dev/shared_gamma_split_20261111/source.sha256).

[The locked input record](../dev/shared_gamma_split_20261111/prepared_v1/INPUTS_LOCKED.json)
links every new discovery/confirmation array,source matrix,proposal lock,
selected identity axis,reference expectation,transformation and fixed seed.

All five original recovery probes and four state controls retain their selected
cells, event genes, reference estimates, masks, coordinates and original count
source. Biological n=0. These already exposed simulations remain regression
mechanism probes. Historical population/feature/nuisance selection remains
conditioned upon and is not independent confirmation.

The full mean and shape are the archived ED+EC and rD+rC after undoing the old
categorical offsets. Exact equality of the old half means/shapes and share1/2
is required. No true simulator dispersion or clone label enters. New proposal
counts are generated once from the original full gene matrix using Philox,
domain20261111/version1,original case seed,global cell ID and original gene
column. Each identity supplies a fixed52-bit midpoint uniform and binomial
inverse CDF. It never evaluates endpoints0/1. All draws and their hashes lock;
there is no redraw, seed search or allocation sweep.

The complete event bins are excluded from new discovery-depth totals. Depth is
centered/scaled on new discovery counts. The initial prevalence uses a median
partition of the new discovery regional count ratio. **Archived q is excluded**:
old discovery observations overlap new confirmation counts. Carrying that q
into a start would violate the new information roles. This adaptation uses no
new confirmation value, hidden truth or favourable outcome.

Preparation saves new discovery/expected arrays separately from new C arrays.
Fitting opens only discovery data and expectations; all24/48/96 fits lock
before C decoding. The frozen null/alternative retain a common linear trend,
Gaussian within-population variance,one shared prevalence,analytic zero-variance
boundary,two/seven starts and the existing optimizer. Fits are successful finite
multistart local solutions, not certified global maxima. No variance floor,
forced copy effect, prior, support reduction or new test threshold is introduced.
Original count conservation D+C=Y must hold through exact global-ID/gene alignment.

The unchanged nine-case continuation criterion requires every true score positive
and numerically resolved with both component masses>=30, and every state score
negative and resolved. The numerical diagnostic remains2×maximum consecutive-
order discrepancy plus rounding budgets. It is not a formal integration enclosure
or confidence interval. A failure prevents full-workflow expansion and any seed,
share or threshold sweep. Passing would merely permit separately frozen full
recovery and fresh complete-workflow null/power, retaining the5% specimen upper
confidence bound and simultaneous power requirements.

## Independent prerequisite

Prepare/numerics23206753 completes0:0 in24seconds on nslab CPU. It locks all
nine data pairs. [Numerical results](../dev/shared_gamma_split_20261111/prepared_v1/NUMERICS.json)
and [independent references](../dev/shared_gamma_split_20261111/checks.py) include:

- Seven exact total-count allocation sums and four exact marginal checks.
  The explicit omitted-tail bound is67/2^66, not a statistical confidence bound.
- A joint2×2 probability minor2/729 and covariance2 for mean4/shape2/share1/2,
  demonstrating positive sampling dependence. No independence claim is inferred
  from RNG replay or empirical zero covariance.
- Exact Fraction conditional population evidence agrees to1.2212453e-15,
  within operation/magnitude budget1.7462298e-10. Treating the new parts as
  independent changes a tested per-cell score by.4071595, far beyond arithmetic
  uncertainty. This is a consequential model mismatch counterexample.
- Exact row/gene permutation,chunk,subset and count-conservation identities.
  Frozen-model rational/adaptive full-PMF gradients,finite differences,nested
  null,component/cell exchange and covariate reparameterization also pass.

These establish the stated computational construction under its assumptions,
not the adequacy of estimated reference/noise parameters or biological validity.
No public package source changed. The521-Python/four-Rust package validation
is historical; a new full package/GPU suite is not claimed.

## Submission checkpoint

All submissions use nslab account/QOS,CPU,one CPU each. Prepare8GB/1h;
nine case jobs8GB/3h;summary4GB/1h;independent reader8GB/1h. Each batch header
was announced before submission. No arrays,interactive QOS,new dependency,
shared installation,running-job change or frozen-output change.

[Prepare command](../dev/shared_gamma_split_20261111/prepare_submission.json),
[nine case commands](../dev/shared_gamma_split_20261111/case_submissions.jsonl),
[summary dependency](../dev/shared_gamma_split_20261111/summary_submission.json),
[reader protocol](../dev/shared_gamma_split_20261111/READER_DESIGN.md),
[independent reader](../dev/shared_gamma_split_20261111/reader.py).

At this checkpoint five cases complete,four are still fitting,and summary/reader
await dependencies. There are no executable failures so far. No whole scientific
gate or recovery improvement follows from partial execution. Complete results,
all attempts and resource accounting will be appended below.

## Complete result:7/9, no promotion

All primary execution completes successfully, but the frozen scientific screen
fails7/9. The [complete summary](../dev/shared_gamma_split_20261111/summary_v1/REPORT.json)
retains all nine cases and all order-specific values.

| Fixed case | Proper shared-noise predictive advantage | Numerical diagnostic | Smaller discovery mass | Decision |
|---|---:|---:|---:|---|
| original20261031 | −.786648 | .00073937 | 191.8693 | Fail |
| original20261107 | +6.142476 | .00079867 | 168.1101 | Pass |
| original20261108 | +2.688681 | .00086132 | 184.1217 | Pass |
| original20261109 | +1.858881 | .00052386 | 208.8343 | Pass |
| original20261110 | +1.471797 | .00069915 | 223.0418 | Pass |
| state202611096001 | +1.47e−7 | .00174490 | 294.4765 | Unresolved;fails screen |
| state202611096002 | −.172322 | .00080104 | 181.4983 | Pass |
| state202611096004 | −.814357 | .00133652 | 128.4366 | Pass |
| state202611096006 | −1.246582 | .00117589 | 292.2123 | Pass |

The original109 contrast changes from the preceding linear-depth score−3.634283
to+1.858881, but original1031 changes from+5.595246 to−.786648. These are
conditional proposal results, not complete clone recoveries. The candidate
retains4/5 true contrasts and3/4 resolved state rejections. The remaining
state result is numerically unresolved, **not** an established false-positive
clone distinction or a successful negative-sign rejection. Every component
exceeds30-cell support; no support rounding can explain either failed criterion.

The same original full RNA counts underlie both experiments, but discovery
allocations, derived covariates, initializers and the joint likelihood differ.
The sign reversal alone cannot isolate a causal benefit of binomial sampling
or establish superiority from one favourable allocation. No redraw, share
change, seed search, selective union of passing cases or threshold adjustment
is justified by this failed experiment. Full-workflow recovery and fresh
null/power expansion are **not run** for this unchanged candidate.

## Independent reader and matched noise diagnostic

[Reader23206791](../dev/shared_gamma_split_20261111/reader_23206791/results/REPORT.json)
reconstructs all27 exact score sums and all9 decisions. All original counts
conserve exactly under global-row/gene alignment. The fixed first/last-ID and
minimum/maximum-depth probes supply36 independent adaptive full-NB-PMF
conditional-score references. Maximum discrepancy5.0598414e-14 passes the
recorded integration/tail/operation budgets. An additional108 binomial draws
are checked using exact integer CDF inequalities,F(d−1)<u<=F(d),with no
floating CDF tolerance. No favourable-score probe selection is used.

The [matched diagnostic](../dev/shared_gamma_split_20261111/ABLATION_DESIGN.md)
then freezes every new count,covariate and fitted parameter and replaces only
joint likelihood assembly with deliberately independent halves. It is an
incorrect joint under the declared shared-noise law, not a competing release
candidate. All three orders and all nine cases remain. Another36 full-PMF
adaptive probes independently verify this diagnostic to5.8730798e-14 maximum
discrepancy. [Complete diagnostic](../dev/shared_gamma_split_20261111/ablation_23206960/results/REPORT.json).

| Fixed case | Proper shared noise | Wrong independent halves |
|---|---:|---:|
| original1031 | −.786648 | −.884076 |
| original107 | +6.142476 | +6.292237 |
| original108 | +2.688681 | +2.412947 |
| original109 | +1.858881 | +2.047301 |
| original110 | +1.471797 | +1.966499 |
| state6001 | +1.47e−7,unresolved | +3.10e−7,unresolved |
| state6002 | −.172322 | −.175103 |
| state6004 | −.814357 | −1.528011 |
| state6006 | −1.246582 | −2.467477 |

No resolved sign or unresolved status changes. Dependence changes numerical
evidence and must be modelled correctly, but this conditional assembly ablation
does not explain the recovery tradeoff by itself. It does not refit parameters,
separate draw effects or estimate a full-workflow false-split rate. The primary
7/9 failure remains unchanged. A mathematically wrong joint would not be
promoted even if a selected accuracy score had improved.

## Complete execution, preservation and limits

All13 jobs use nslab account/QOS,CPU,one CPU. No running job is disturbed.
Prepare8GB/1h,cases8GB/3h,summary4GB/1h,reader8GB/1h,diagnostic8GB/1h.
Every submission batch header was announced before execution.

| Job / role | Progress | State | Elapsed | Peak process RSS,KiB |
|---|---|---|---:|---:|
| 23206753 numerics/preparation | 1/1 ████████ | Complete0:0 | 24s | 206108 |
| 23206763 original1031 | 3/3 ████████ | Complete0:0 | 1m47s | 272896 |
| 23206764 original107 | 3/3 ████████ | Complete0:0 | 4m11s | 260488 |
| 23206765 original108 | 3/3 ████████ | Complete0:0 | 1m40s | 262996 |
| 23206766 original109 | 3/3 ████████ | Complete0:0 | 1m17s | 256160 |
| 23206767 original110 | 3/3 ████████ | Complete0:0 | 1m03s | 261752 |
| 23206768 state6001 | 3/3 ████████ | Complete0:0 | 5m23s | 289456 |
| 23206769 state6002 | 3/3 ████████ | Complete0:0 | 10m53s | 294152 |
| 23206770 state6004 | 3/3 ████████ | Complete0:0 | 15m23s | 291604 |
| 23206771 state6006 | 3/3 ████████ | Complete0:0 | 1m36s | 294780 |
| 23206772 complete summary | 9/9 ████████ | Complete0:0 | 5s | 41504 |
| 23206791 independent reader | 36/36 ████████ | Complete0:0 | 12s | 113860 |
| 23206960 matched diagnostic | 9/9 ████████ | Complete0:0 | 42s | 526628 |

Running0,pending0,complete13,executable failures0. Scientific failures2/9,
including one unresolved state control. Peak process RSS526628KiB is measured
by GNU-time; it is not the allocation or aggregate concurrent memory. The
longest fit remained on its original job and completed; no observation timeout
or slow fit triggered a restart, duplicate or changed draws.

[Accounting](../dev/shared_gamma_split_20261111/accounting_complete.tsv),
[process resources](../dev/shared_gamma_split_20261111/resources_complete.json),
[reader command](../dev/shared_gamma_split_20261111/reader_submission.json),
[diagnostic command](../dev/shared_gamma_split_20261111/ablation_submission.json).

Every executed source and reference module retains its hash. The copied
linear-depth family and numerical checks are unchanged. All91 previously
validated package/native/data/test identities and the earlier arm-scoring
correction remain exact. [Integrity](../dev/shared_gamma_split_20261111/source_integrity_complete.json).
No shared environment,dependency,package,default,existing test,scientific cutoff
or frozen result changed. The ledger adds this complete failed experiment and
retains all282 earlier entries. The earlier partial checkpoint remains above
as history; it is superseded by the complete accounting here.

Computational correctness passes the stated tests. Statistical identification
remains conditional on working reference/noise/state laws and historic selection.
Biological validation remains absent. The unresolved state control must not be
relabelled a measured false-split event. Original production recovery is still
2/5,continuous-state false splits14/200,and correlated selected-interval
coverage54/1000. No cell-membership accuracy,spatial profile,bulk accuracy,
competitor superiority or biological5% error claim is established here.
The full scientific/SOTA objective remains active and unachieved.
