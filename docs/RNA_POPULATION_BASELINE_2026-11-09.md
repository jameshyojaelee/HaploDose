# Fitted RNA populations do not solve the clone error/power problem

Follow-up: the original67 literal MAP-change cases include4 complete label
reversals, leaving63 changed partitions. None affects an accepted evaluation
distinction. The largest score change is an improvement, not a loss of a coarse
solution. [Corrected independent audit](RNA_POPULATION_OPTIMIZER_2026-11-09.md)
supersedes the numerical interpretation below; all original results remain.

A simple genome-wide RNA population comparator fails the unchanged error/power
goal. It accepts8/200 matched-NB null distinctions and10/200 state-null
distinctions, with simultaneous upper bounds11.02% and12.47%. It detects14/200
two-population alternatives. Do not promote it or describe fitted RNA mixtures
as validated genomic clones. Original recovery2/5 and all frozen failures remain.

## Scientific contract and the removed oracle inputs

The estimand is a distinguishable population in observed, normalized regional
RNA expression, plus assignments conditional on that working model. The
comparator normalizes observed counts by their observed library totals, forms
half-count log ratios to the same known reference supplied to HaploDose,
centres features across cells, and fits a tied-variance two-Gaussian mixture to
the leading full-SVD component across all240 bins. A single fitted Gaussian is
the competing RNA-state model. This Gaussian model on transformed RNA is a
working approximation, not an NB count likelihood or biological DNA model.

Unlike the earlier count oracle, fitting receives no true depth, event region,
CNA effect, state distribution, simulated label or balanced class size. The
generating reference remains known and shared with HaploDose: this is a
conditional known-reference experiment, not a test of reference identifiability.
The machine-readable `oracle_inputs_used=false` is limited to those excluded
membership-oracle inputs; it does not establish that reference construction is
solved in biological specimens.

Genomic clone interpretation would require independent same-cell DNA, or
independently justified restrictions excluding regional transcriptional-state
populations. Mixture components alone cannot distinguish those biological
explanations. There are no biological specimens or new DNA measurements.

## Locked complete comparator workflow

Reuse all998 exposed observations from the paired joint experiment, requiring
exact count, reference, dispersion, truth, state, mean and reset hashes. Each
condition has600 cells and240 bins. There are199 calibration specimens for each
of the matched-NB and correlated-state nulls, and200 evaluation specimens in
each null and the two-CNA alternative. These are simulated specimens, not
independent patients or new evaluation cohorts.

PCA, nuisance fitting, mixture selection and MAP support checks are repeated
in every null and observed workflow. The statistic is twice the mixture
log-likelihood improvement over a single Gaussian. Both groups must meet the
existing30-cell MAP minimum, applied in both observed and calibration scores.
Use the maximum of the two inclusive null-bank rank p-values at the single
prespecified0.05 operating point. No chi-squared approximation, BIC selection,
gate search or post-outcome resampling is used. All original criteria remain.

Independent Fraction ranks agree with strict order-statistic decisions, including
ties. The new analysis includes156 simultaneous one-sided95% exact bounds:
the previous153 endpoints are recomputed in a separate artifact alongside the
three new endpoints. Historical153-family results are unchanged. Independent
80-digit binomial inversion and SciPy beta references agree exactly for the
two null bounds and within6.94e-18 for power. Bounds condition on the realized
calibration banks and the stated simulation laws.

| Evaluation condition | Accepted distinctions | Rate | Simultaneous bound |
|---|---:|---:|---:|
| Matched NB, one CNA population |8/200|4.0%|Upper11.0230%|
| Correlated state, one CNA population |10/200|5.0%|Upper12.4699%|
| Two CNA populations |14/200|7.0%|Lower2.3559%|

Both null upper-bound requirements fail. This is not itself proof that either
true rate exceeds5%. Nominal rank control is not a passed conditional empirical
upper-bound criterion. Low power also prevents qualifying a mostly unsplit
procedure as reliable clone recovery. No extra draws or gate adjustment follows.

## Membership and numerical uncertainty

Before the specimen gate, two-group all-cell correctness is60.3508%, using only
postfit label-permutation matching. Accepted-case correctness is65.4643%, but
accepted two-group coverage is only7%; correctness over all cells retaining
abstention is4.5825%. All600 cells remain in each denominator. This baseline
hard-assigns cells before gating and has no normal/reference-like classification
or calibrated individual membership confidence. Its detection count is always
an exactly-two result by construction; that restriction is not proof of accurate
membership. Do not substitute these metrics for the original recovery assertion.

The earlier67.26% count reference uses unavailable generating information. The
new60.35% result does not close that gap or isolate how much is due to projection,
nuisance estimation, model approximation or finite optimization. Comparisons to
HaploDose's fragmented partitions remain descriptive, not superiority tests.

The original EM solver passed its small preflight but exhausted2000 iterations
in both full runs. Preserve those failed jobs and outputs. The repair directly
optimizes the same likelihood using analytic gradients checked by independent
finite differences. Component means are bounded by the observed range; shared
variance is bounded by the optimal hard two-means residual variance and total
variance. These follow from likelihood-optimum moment identities, rather than
an arbitrary covariance ridge. Ten deterministic starts and the iteration budget
remain fixed. A converged selected solution is required; global optimality is
not certified. Independent scalar SciPy normal-density sums check every fit.
Canonical count-row order makes score, loadings and memberships exactly invariant
to a tested cell permutation.

Refinement from optimization targets1e-8 to1e-10 changes no600 evaluation rank
decisions in any coarse/fine observed-bank combination. However67 of998 cases
change MAP assignments, and the maximum score change is6.59098. These are
nontrivial numerical sensitivity findings, not a convergence certificate. Do
not present memberships as numerically settled or promote the solver. A useful
follow-up would retain the better coarse solution during continuation and check
independent optimization on implicated cases before any full corrected replay.
That numerical work cannot supply missing biological identifying information.

## Execution, preservation and decision

| Job | Stage | Seconds | Process peak RSS KiB | Exit |
|---|---|---:|---:|---:|
|23062501|original EM preflight|21|213908|0|
|23062502|original calibration|29|226240|1|
|23062511|original evaluation|24|229524|1|
|23063191|direct-optimizer preflight|12|203736|0|
|23063192|398 calibration fits|266|217836|0|
|23063193|600 evaluation fits|396|220180|0|
|23063368|initial three-endpoint summary|22|199716|0|
|23063621|full156-family summary/audit|23|200660|0|

Final evidence is `dev/rna_population_baseline_20261109/summary_v2/report.json`;
the prior summary lacks the retained153 recalculations and is superseded for
family-wide reporting. Designs precede successful execution. Immutable source,
per-case hashes/probabilities/loadings/optimizer diagnostics, environment/backend
records, failed attempts and resources remain under the same task directory.
Scripts: `benchmarks/hardening/rna_population_baseline.py` and
`benchmarks/hardening/rna_population_summary.py`.

All jobs execute on Intelv4 CPU nodes with innovation/nslab; no shared installation
or new dependency. Static hash audit verifies73 package/native/test/README files
unchanged from the turn baseline. The full package suite was not rerun because
production is unchanged; the executable comparator and audit were run on compute.
Do not promote this method, dilute the5% bound, reinterpret RNA mixtures as DNA
clones or claim SOTA. The missing biological null panel remains n=0, and the
original recovery, uncertainty and false-split/power failures remain open.
