# Known-nuisance cell-population prediction

The scientific goal remains useful genomic profiles and supported clone
distinctions from RNA, with evidence at the claimed resolution. This experiment
addresses the review's findings 2, 4, 5, 7 and 9: whether count nuisance and the
prediction target jointly cause the exposed true-contrast rejection. It does
not establish a genomic interpretation or achieve the specimen error target.

## Frozen comparison

Use the same exposed seed 202611099001, 600 cells, native selected chr12 event,
raw gene counts, original gene vetoes and deterministic cell roles as the
[completed non-oracle cell-population test](CELL_POPULATION_PREDICTION_2026-11-11.md).
Supply the already qualified known baseline, NB size 10 and full affected-region
compositional share from the [count-nuisance diagnostic](STATE_NUISANCE_ORACLE_2026-11-11.md).
Keep depth covariate zero because exposure is supplied. No truth labels or
generating amplitudes enter fitting. This is the last missing combination of
these two fixed mechanism interventions, not a search over new models.

| Prediction target | Estimated count nuisance, state / true score | Known count nuisance, state / true score |
|---|---:|---:|
| Further counts from the same cells | −1.950379 / −2.457458 | −3.062297 / −2.270844 |
| Other cells, fixed 300/300 roles | −0.107550 / −2.750987 | −0.188559 / −5.645948 |

The scores have different conditioning sets; comparing their numerical
differences as factorial main effects would be incorrect. The prespecified
paired success criterion is common: a resolved negative state score, and a
resolved positive true-contrast score with both original mass/MAP minimums of
30, numerical membership resolution and both original gene vetoes. No choosing
the favourable role, seed, threshold or older criterion.

The working observation model is shared-Gamma gene NB counts with Gaussian
regional state and optional two-population shift. The estimand is predictive
support for a distinct regional RNA population, not DNA identity. Competing
regional transcription can generate the same observations. A pass identifies
joint requirements of this oracle mechanism; a failure falsifies their
sufficiency. Neither outcome estimates population false-split rate or power.
Biological n=0; independent base seeds=1.

## Implementation and independent evidence

[Frozen design](../dev/known_cell_population_20261111/DESIGN.md),
[executed study source](../dev/known_cell_population_20261111/study.py),
[independent reader](../dev/known_cell_population_20261111/reader.py),
[source identities](../dev/known_cell_population_20261111/SNAPSHOT.json).

The reader verifies both prediction locks before opening truth; reconstructs
raw counts and gene/cell correspondence, known baseline and region denominator;
independently checks complete-PMF training/target memberships, adaptive
integration at consequential cells, exact Fraction score sums and exact gene
votes. Original 24/48/96 quadrature, model fits, starts, support, numerical
diagnostics and scientific criteria remain. Shared selection and oracle
nuisance prohibit independent-confirmation or deployable-calibration claims.

## Execution

| Job | Purpose | Allocation | Dependency |
|---|---|---|---|
| 23213092 | State condition | CPU Intel-v4, 1 CPU, 8 GB, 3 h | None |
| 23213093 | True contrast | CPU Intel-v4, 1 CPU, 8 GB, 3 h | None |
| 23213094 | Independent reader | CPU Intel-v4, 1 CPU, 8 GB, 1 h | afterok both fits |

All account/QOS nslab. Existing isolated Conda environment, neutral execution
directories, source hashes before/after, time-v resource records; no arrays,
new dependencies or shared installations. Before submission: AST and bash
syntax passed, all 91 protected identities verified, immutable source copied.
The original dirty tree, frozen results, prior claims and running/unrelated
jobs remain unchanged.
[Exact submission commands](../dev/known_cell_population_20261111/submissions.json).

Results and executed validation are pending at this checkpoint. No public
default or readiness claim changes.

## Completed result: combined intervention is insufficient

All three jobs complete 0:0. The independent reader verifies the failed paired
criterion. State score is −0.18855881254170187 with numerical diagnostic
0.00037297960247238064; true-contrast score is −5.645947753250283 with diagnostic
0.0009560615062677218. Both reject. This is not a near-boundary numerical outcome.

| Condition | Two directional component masses | Independent minimum support lower bounds | Original gene votes |
|---|---|---|---|
| State | 148.839/151.161; 32.758/267.242 | 148.839252; 32.757987 | 58/92 and 64/91, both pass |
| True contrast | 202.834/97.166; 42.208/257.792 | 97.165767; 42.208194 | 81/93 and 76/93, both pass |

No mass floor or gene veto explains the true-contrast rejection. Both component
orientations and all numerical memberships resolve. The unchanged nominal
gene-sign test also passes the continuous state; it is not genomic confirmation.
Raw classifiers score 342/600 and 343/600 respectively, with zero numerical
abstentions, but are not applied because the selection score rejects.

Independent validation covers 1200 complete-PMF training/target probabilities
per condition, all 18 exact rational score sums, 12 adaptive full-PMF probes
(largest discrepancy 1.3322676295501878e−14), exact raw-count and gene/cell maps,
known baseline/size/region denominator, role identities and every locked output
hash. Numerical checks pass without a tolerance change. This establishes correct
computation of the negative result; it does not establish a global optimisation
optimum, finite-count information limit or population error rate.

Known nuisance and prediction-target changes are insufficient separately and
together under this working procedure. Reject their promotion as the repair.
[Inference] Remaining possibilities include finite-count signal, estimation
cost and the fitted population family/optimisation. Their relative contributions
are not identified by this one paired draw. The existing noiseless-amplitude
oracle cannot substitute for a finite-count experiment. No new calibration bank
or parameter/role search is justified by this negative result.

Fit runtimes are 8m16s and 8m06s; the reader completes in 21 seconds by scheduler
accounting. Scheduler fit peaks are 237276/237392 KiB; time-v process records
are retained separately. All account/QOS nslab; no failures or retired jobs.

[Authoritative independent report](../dev/known_cell_population_20261111/reader_23213094/results/REPORT.json),
[state resource record](../dev/known_cell_population_20261111/job_23213092/resources.txt),
[true resource record](../dev/known_cell_population_20261111/job_23213093/resources.txt),
[reader resource record](../dev/known_cell_population_20261111/reader_23213094/resources.txt).

No public/default changes, no weakening of the false-split/power or recovery
criteria, and no change to frozen outputs. Biological n=0. Research recovery
remains 5/5, recorded production 2/5, state errors 14/200 and selected correlated
coverage 54/1000. The full scientific goal remains unachieved.
