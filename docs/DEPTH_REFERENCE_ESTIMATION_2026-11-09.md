# Estimated reference depth effects: benefit depends on support and transport

The conditional intercept correction reduces RNA-profile error with well-supported
reference depth variation, but is unsafe as a universal correction. Estimated
slopes amplify noise under narrow support and introduce effects when the true
slopes are zero. Even known slopes cannot transport a reference-specific
expression programme to a query lacking that programme. No production correction
or calibration claim is promoted. Frozen evidence and criteria are unchanged.

## Scientific question and experimental contract

The estimand is pooled expected normalized RNA abundance at a specified query
exposure distribution, compared with same-genotype references. It is not DNA
copy number. The observation model generates raw integer binned RNA counts using
independent gamma-Poisson draws with variance mean+.05*mean². Exposures and the
depth-dependent compositions are known generating inputs. The fitted slope uses
observed count totals, which differ from those exposures; the query prediction
uses the known exposures in this conditional diagnostic. This is not a deployed
end-to-end uncertainty experiment.

Reference depth association must represent the same nuisance mechanism in the
query. The strongest competing explanation is a depth-correlated reference
expression programme. The confounded challenge changes that programme between
reference and query while preserving genotype. Heldout references share the fit
reference programme, so successful heldout prediction cannot establish transport.
Independent corresponding RNA/DNA and depth-matched biological replication would
be needed to distinguish a genomic explanation in real material.

Eight independently seeded simulation replicates in each of six prespecified
conditions are conditional computational replication, **biological n=0**. Fit200,
heldout100 and query100 cells per scenario, except30 fit references in the small
condition. Forty bins; fixed known RNA gain0.6*log1.5 in the first10. The programme,
depth ranges, sample sizes, seeds and comparison family were frozen in
[the count design](../dev/depth_reference_estimation_20261109/DESIGN.md).
No thresholds or cases were selected using outcomes.

## Count-panel results

Mean null centred log-rate MSE over all eight replicates:

| Reference/query condition | Current expectation | Estimated-slope intercept correction | Known-reference-slope correction |
|---|---:|---:|---:|
| Adequate depth range, query within range |0.00399348|0.00026261|0.00025641|
| Adequate depth range, query beyond range |0.30930058|0.00639593|0.00029480|
| Narrow depth range, query beyond range |0.38793159|0.62396112|0.00024891|
| Only30 fit references |0.33021826|0.33021826|0.00202382|
| True slopes zero, query beyond range |0.00033397|0.00580445|0.00030965|
| Reference-specific programme |0.33685558|0.18840044|0.18757702|

The known-slope comparator is an oracle mechanism check, not a deployable method.
Under the supported challenge, mean first-region gain is0.23158 after estimated
correction versus generating0.24328. Under extrapolation it is0.20662. Narrow
support reverses it to−0.39718; confounding reverses it to−0.61144. Thus a general
gain-preservation claim is rejected. These are RNA contrasts, not clone recovery
rates, DNA sensitivity or biological performance.

The production slope routine returns all zeros with30 references because its
default minimum30 requires at least60 cells for quintiles. This is lack of
estimating support, not evidence of no depth effect. Under adequate reference
support, fitted noisy slopes still create appreciable extrapolation error when
the true slopes are zero. Intercept correction alone does not fix that problem.

## Computational evidence

[The panel](../benchmarks/hardening/depth_reference_estimation.py) uses the isolated
validated wheel and actual installed slope/expectation routines. It checks
positive/conserved query expectations, reference reproduction, coarse/fine
calibration, reference permutation and a heldout-only mutation. The mutation
leaves estimated slopes and query predictions exactly unchanged. Query targets
are constructed after each predictor is written; oracle values are never used
to choose fitted settings.

[The independent audit](../benchmarks/hardening/depth_reference_estimation_audit.py)
reconstructs quintile boundaries through explicit order statistics, pools counts
with scalar arithmetic, recomputes softmax means and reference equations, and
checks all384 contrasts/1152 scalar metric values. Operation-count precision
bounds apply to numerical equality only; they are not scientific error criteria.
All48 cases pass. Full tables and source/environment identities are retained in
[the task directory](../dev/depth_reference_estimation_20261109).

The first launch failed because this diagnostic passed sparse counts to an
internal dense-count interface. Its output and failed exit are retained; the
replacement changes only that call. No production code or original tests changed.

## Reference uncertainty follow-up

The separately frozen
[uncertainty design](../dev/depth_reference_estimation_20261109/UNCERTAINTY_DESIGN.md)
resamples fit references199 times per scenario, refitting the full nuisance
estimator. It compares a simultaneous whole-profile band and a fixed first-region
interval with the known targets and gain. These conditional RNA events are not
the5% whole-specimen clone endpoint. All9552 bootstrap replicates finish without
failing the reference-residual checks; this does not assert all solver status
flags indicate convergence or that their statistical model is correct.
The two intervals target different functionals: a maximum-bin simultaneous band
and the mean of the first10 separately median-centred log rates. Those contrasts
are not substituted for the earlier median-centred ratio estimand.

| Condition | Whole-profile null misses /8 | Fixed-region null exclusions /8 | Gain positive lower bound /8 |
|---|---:|---:|---:|
| Supported |0|1|8|
| Extrapolated |0|0|8|
| Narrow |2|6|0|
| Small |8|8|0|
| Zero slopes |0|0|8|
| Reference-specific programme |8|8|0|

Sampling intervals do not repair absent support or nontransportable biology.
The small-reference slope remains zero in every resample; its narrow reference
sampling interval cannot represent the missing depth effect. Under narrow
support, whole-profile bands are large (median2.203 log units) yet still miss2/8;
the fixed-region interval misses6/8. These observations falsify a general
uncertainty-repair claim and show why both coverage and retained signal matter.
For zero slopes, wider intervals expose estimated-slope noise in these eight
cases; this is not evidence that introducing that noise is beneficial.

Even0/8 gives an exact one-sided95% upper null-miss bound0.31234, or0.45033
with Bonferroni simultaneous bounds over the six conditions. Two/8 gives0.59969
and0.71660. These bounds describe conditional simulations, not biological
specimens or clone error. The nominal95% bootstrap quantile is an approximate
resampling procedure with199 draws, not an exact coverage guarantee.

[Independent interval audit](../benchmarks/hardening/depth_reference_uncertainty_audit.py)
recounts all96 quantiles,288 decisions,12 exact-binomial bounds by polynomial
inversion, and144 heldout deviance scores by scalar arithmetic. All pass. The
heldout score is a Poisson composite predictive diagnostic for these NB-generated
counts; it is not claimed to be their joint likelihood.

## Compute and preservation record

| Job | Purpose | Account | Exit | Runtime |
|---|---|---|---:|---|
|23081016|Initial diagnostic, sparse input failure|innovation|1|15s|
|23081037|Count panel with corrected interface|innovation|0|19s|
|23081229|Independent count-panel audit|nslab|0|3s|
|23081679|9552 reference-resampling fits|innovation|0|3m12s|
|23081910|Dependent interval and score audit|nslab|0|2s|

All ran on allocated CPU nodes from immutable snapshots, using the previously
validated noneditable wheel. Detailed resources, hashes, dependencies, source
identities and retained failure outputs are in
[jobs.json](../dev/depth_reference_estimation_20261109/jobs.json) and linked run
manifests. The original package/test/native/README files remain unchanged.

## Decisions and next implementation gate

Computationally, the intercept correction is reproducible and addresses the
known-origin example. Statistically, unknown slopes, insufficient depth support
and nontransportable reference programmes remain limiting assumptions.
Biologically, no DNA claim is validated by these simulations.

Do not implement an expected()-only patch. A consistent implementation must
carry reference depths/origins through reference construction and mixture fitting,
initial pooled fitting, per-cell rates, profile expectations and confirmation
simulation. Their intercept and pooled-rate estimands differ when depths vary.
Reference uncertainty and an explicit unsupported/abstention path must be tested
before default promotion. A heldout score alone cannot certify transportability.

The original three recovery failures, complete-workflow false-split/power failure
and uncertainty failures remain open. No result here changes those verdicts,
admits a new evaluation cohort, or establishes superiority over simple baselines.
