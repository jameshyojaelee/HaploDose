# Reference averaging improves bulk point accuracy, with a loss-detection cost

[Observed pattern, not a guarantee] The implemented fixed reference ensemble
raises mean bulk profile Pearson to0.770411 from the retained0.764817, versus
0.732762 for equally ensembled pseudobulk. It leads all84 methods in this exposed
regression panel and preserves all45 preceding score tables exactly. Loss recall
falls0.602747→0.540476 and FPR rises0.043737→0.044760. The fixed conditional
promotion gate fails; retain a research profile candidate, not a default or
independently validated superiority claim. All3 CPU jobs complete0:0.

## The implemented improvement and its scientific contract

Reference perturbations already changed the preceding model's performance.
Average their predictions inside the estimator rather than merely reporting
reference sensitivity afterward. Use ALL8 existing frozen reference-bootstrap
fits001..008, excluding central fit000, with an equal arithmetic mean. Apply the
IDENTICAL operator to every39 method with those resamples, including pseudobulk.
No reference, draw, number of draws, mean/median operator or patient is selected
by DNA. The primary is reference_bag_averaged_gaussian_reference_median, fixed
before any ensemble score. [Frozen design](../dev/reference_bagged_profile_20261110/DESIGN.md).

Estimand: regional relative RNA abundance at the fixed validation-window
resolution, averaged under these empirical reference perturbations. This
finite8-fit ensemble is a deterministic estimator, not a converged posterior
integral. Observation model remains the existing RNA/Gaussian regional models
inside each fit. Reference sampling must represent a relevant reference
population for a variance-reduction interpretation. Systematic lineage
composition or regional transcription can persist in every member and explain
DNA-shaped profiles. Independently corresponding DNA at matched resolution
tests genomic prediction; repeated independent RNA references distinguish
reference robustness from averaging a shared bias. This former evaluation
cohort stays regression-only. New biological replication=0; no DNA copy-state,
clone membership or spatial validation follows.

The smallest decisive test is identical aggregation of candidate and baselines
on the existing8 fits, then unchanged patient-level scoring. Falsification:
the fixed primary must beat ALL83 alternatives in paired mean Pearson before
an uncertainty continuation; no cherry-picking another ensemble. The RNA stage
reads only method arrays and sample IDs. Source archives also contain frozen DNA
arrays, which that stage does not access. All39 outputs/hash records lock
before the dependent scorer reads DNA. Exact Fraction mean and permutation
checks pass; all8 finite-support masks and specimen orders must agree.

## Complete matched performance

All45 preceding methods, including the failed variance-composition ablations,
remain. Add39 identically aggregated methods for84 total. Keep all48 patients,
24/24 calibration/evaluation split,23 informative correlations and13 recall
patients, masks, scorer, calibration threshold rule and empty-output handling.
Every45 earlier patient-score table reproduces EXACTLY.

| Method | Mean Pearson /23 | Loss recall /13 | Neutral loss FPR /23 | Precision /calling patients |
|---|---:|---:|---:|---:|
|Primary reference ensemble|0.7704112303|0.5404762|0.0447602|0.4240196 /17|
|Previous additive boundary average|0.7648167082|0.6027473|0.0437372|0.4837963 /18|
|Equally ensembled strong pseudobulk|0.7327616523|0.3996795|0.0419880|0.3235556 /15|
|Original strong pseudobulk|0.7302613671|0.3996795|0.0412445|0.3585714 /14|

Primary paired mean Pearson advantages:0.00559452 over the previous best
(15wins/8losses),0.03764958 over equally ensembled pseudobulk(22/1),
0.04014986 over original pseudobulk(22/1). Its closest point alternative is
ensembled hard Gaussian reference median0.7678296, gap0.00258161.
All83 point contrasts are positive, so the prespecified point gate passes.
[Full84-method report](../dev/reference_bagged_profile_20261110/score_v1/report.json),
[all patient scores](../dev/reference_bagged_profile_20261110/score_v1/patient_scores.tsv).

The recall/FPR tradeoff prevents describing this as a general loss-detection
improvement. Precision populations differ and remain caller-conditional;
their marginal means are not a paired precision result. The candidate currently
exists at window resolution in the research implementation, not a newly
delivered public bulk API. A full-resolution export is not claimed.

## Conditional uncertainty and the failed continuation gate

After point results, freeze the conditional uncertainty procedure BEFORE its
draws. Resample all24 evaluation patients jointly20000 times, seed2026111031;
retain paired missingness. Test primary against ALL83 alternatives for Pearson,
loss recall and loss FPR(lower favorable):249 fixed contrasts. No comparator
or unfavorable endpoint disappears. Same paired mean is effect, test and
interval estimand. Use basic bootstrap Bonferroni intervals, centered-error
two-sided p, Holm over249, and exact-binomial simultaneous Monte Carlo brackets.
Independent rational interpolation, Holm and actual-draw fsum checks pass;
maximum arithmetic discrepancy5.55e-17.
[Conditional design](../dev/reference_bagged_profile_20261110/CONDITIONAL_DESIGN.md).

| Comparison | Paired mean Pearson gain | Approximate simultaneous interval | Upper-MC Holm p |
|---|---:|---:|---:|
|Previous best|0.00559452|[-0.00211158,0.01222393]|1.0|
|Equally ensembled pseudobulk|0.03764958|[0.01049451,0.06431357]|0.128264|
|Ensembled hard Gaussian median|0.00258161|[-0.00116320,0.00557832]|1.0|

These approximate interval point estimates do NOT include their Monte Carlo
brackets. At the249-family tail,20000 draws cannot resolve many extreme
quantiles:72 comparisons remain numerically unresolved, including Pearson
against pseudobulk. All lower-endpoint conservative MC brackets extend to
negative infinity. Retain those limitations rather than reporting the positive
point interval as a passed criterion. Zero contrasts pass the combined
interval/upper-MC Holm rule. Profile and joint gates both fail.

The failure against the previous best and close ensembled hard model is
decisive even using the favorable end of their lower-endpoint MC brackets
(-0.00164797 and-0.000654898). Thus additional draws to resolve pseudobulk alone
cannot rescue the ALL83 gate. No larger Monte Carlo run is selected to produce
a winning statement. Loss recall versus the previous best also has negative
paired point effect-0.0622711. Stop expensive outer-reference/inner-ensemble
uncertainty expansion; preserve the point gain and all negative endpoints.
[All249 effects, intervals, MC brackets and decisions](../dev/reference_bagged_profile_20261110/conditional_v1/report.json).

This bootstrap is conditional on the frozen8-fit estimator and calibration
thresholds. It does NOT include full reference-estimator or calibration
uncertainty, establish coverage, or replace the earlier failed69-family joint
analysis. Neither a conditional pass nor additional compute would establish
independent biological validity. The statistical result does not prove no
population benefit; it withholds the specified superiority claim.

## Every execution, code and preservation

Research sources:
[reference_bagged_profile.py](../benchmarks/hardening/reference_bagged_profile.py),
[reference_bagged_uncertainty.py](../benchmarks/hardening/reference_bagged_uncertainty.py).
Two immutable snapshots, read-only validated wheel, neutral directories,
protected frozen paths, all input/source/environment/backend/native identities
and resources are retained. No original test, production/native code, dependency,
default, criterion, reference fit or frozen artifact changes.

| Job | Role | Account/QOS | CPU / RAM / requested time | Actual time |
|---|---|---|---|---|
|23157259|all39 RNA ensembles|innovation|1 /4GB /1h|7s|
|23157260|dependent84-method scoring|nslab|1 /4GB /1h|6s|
|23157280|249-family conditional uncertainty|nslab|1 /8GB /1h|12s|

All CPU/v6&amd,0:0, no arrays/interactiveQOS. The short tasks have scheduler RSS
recorded as0; actual `/usr/bin/time -v` peak RSS is286564KB and complete
resources remain in analysis job records.
[Exact commands](../dev/reference_bagged_profile_20261110/jobs.json),
[accounting](../dev/reference_bagged_profile_20261110/accounting_final.tsv),
[resources](../dev/reference_bagged_profile_20261110/resource_summary_final.json),
[preservation](../dev/reference_bagged_profile_20261110/preservation_final.json).
No new full package/GPU suite was run for research-only changes. Production
clone recovery stays2/5; the prior5/5 research recovery still fails state-null
error control. The5% whole-specimen confidence/power target, all three original
recovery failures and full scientific goal remain unresolved.
