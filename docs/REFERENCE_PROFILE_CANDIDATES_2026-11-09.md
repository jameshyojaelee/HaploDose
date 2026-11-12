# Reference weighting and adaptive RNA profiles: no strongest-baseline win

Three implemented profile candidates were tested on the two locked Kinker
mechanism cases. The simplest continuous partition improved both cases, but
its frozen181-line regression replay did not beat the strong smoothed pseudobulk
baseline consistently:75 wins/106 losses, mean paired Pearson change+.000406,
median-.006229. Reject a default replacement. This is a new matched analysis,
not an overwrite of historical results or independent evaluation.

## Measurements, estimand and fixed comparisons

The estimand is a regional RNA log-abundance contrast relative to the declared
fit-reference population. Independently measured relative DepMap DNA is used for
prediction assessment at1Mb support; it does not validate cell membership or
absolute DNA dosage. Kinker integers were reconstructed from published CPM and
are not established raw molecule counts. No count-calibration claim is made.

Exact query/ref IDs, selected genes, mappings and bins come from the locked
[primary measurement-chain record](../dev/kinker_measurement_chain_20261109/chain_v1/identities.json).
The two cases are2313287_STOMACH(rotation0,103 query cells/2118 fit references)
and42MGBA_CENTRAL_NERVOUS_SYSTEM(rotation1,517/2467). They were selected previously,
not from these candidate outcomes. The inputs use10133/10116 RNA-selected genes;
DNA support retains10019/10003 genes and2027/2024 common1Mb windows.

All arms use median-centered log((query+.5)/(expected+.5)), where expected is
pooled reference share times the query retained total. The pseudocount is the
existing comparator's rule, matched across arms. Normalize over all RNA-selected
genes before restricting DNA support. Reference DNA is the query line minus
UMI-weighted fit-reference line DNA. RNA and DNA cultures differ, and published
DNA values remain uncertain reference measurements.

Common genomic aggregation, genes and smoothing are mandatory. The existing
five-bin running median resets at hg38 chromosome arms and any declared genomic
blocks in both arms. Historical bins lack explicit gap blocks; this does not
claim newly resolved physical gaps. New normalization/support/reset code differs
from the old181-line scorer, so do not compare the new.668-scale correlations
directly with the historical.708 summary. That discrepancy has not been isolated
to a single cause. Within this experiment all candidate/baseline comparisons
share the same scorer and support. CopyKAT/infercnvpy were not rerun; no new
competitor-family success verdict is issued.

## Implementations and successive fixed designs

[Research implementation](../benchmarks/hardening/reference_weighted_pseudobulk.py)
commits and hashes RNA predictions before loading DNA values.

1. **Total reference variability:** pool each declared reference cell line,
   calculate per-gene log shares with the same pseudocount, estimate across-line
   sample variance(ddof1), and use inverse variance in each existing RNA bin's
   weighted median. Nonfinite/nonpositive variance fails; no fitted floor.
   Exact half-weight ties average adjacent values, reproducing ordinary even-n
   medians at uniform weights. The assumption is that reference-variable genes
   are less reliable dosage reporters. Reference CNAs and lineage effects are
   competing explanations for that variability. Corresponding independent DNA
   predictions distinguish practical benefit, rather than the weight formula.
2. **Within-bin variability:** subtract each reference gene's leave-one-gene-out
   local median before calculating variance. This makes residuals invariant to
   a common per-line regional shift. Singleton bins receive unit weight because
   their one-gene median cannot depend on weight; it is not a variance estimate.
   An exact-arithmetic fixture checks shift invariance. This fixes the proposed
   weighting mechanism's invariance, not its biological interpretation.
3. **Continuous partition:** reuse the existing exact weighted MDL partition
   with SSE+3log(n)*(segments-1) per genomic block. One log(n) penalizes a new
   continuous mean, two log(n) an unspecified breakpoint. Estimate a constant
   working SD from within-block adjacent-bin differences as
   MAD/Phi^-1(.75)/sqrt(2), withPhi^-1(.75)=.6744897501960817. Require a positive
   finite scale; no floor or DNA-selected penalty grid. Sparse breakpoints and
   locally stationary Gaussian residuals are assumptions. This is RNA smoothing,
   not a calibrated CNA discovery or interval-coverage procedure. Correlated
   expression programmes can also form smooth regions.

[Initial design](../dev/reference_weighted_pseudobulk_20261109/DESIGN.md),
[within-bin follow-up](../dev/reference_weighted_pseudobulk_20261109/WITHIN_BIN_DESIGN.md),
[partition follow-up](../dev/reference_weighted_pseudobulk_20261109/PARTITION_DESIGN.md)
and [frozen regression design](../dev/reference_weighted_pseudobulk_20261109/REGRESSION_DESIGN.md)
record when each mechanism was tested. Each retains both cases and the failed
preceding candidate; there is no per-DNA-case winner selection.

## Locked-case results, identical support and operators

| Method |2313287 Pearson|42MGBA Pearson|
|---|---:|---:|
| Uniform pseudobulk bin medians |.625235|.715421|
| Uniform +sm5 |.748763|.776089|
| Total-variance weights |.591462|.742251|
| Total-variance weights +sm5 |.724526|.801313|
| Within-bin weights |.606127|.733609|
| Within-bin weights +sm5 |.721928|.798735|
| Uniform +MDL |.766858|.810360|
| Within-bin weights +MDL |.726216|.837619|

Neither weighting rule improves both cases under matched smoothing. Uniform MDL
does, so freeze that simpler candidate and omit weighting from the broader
replay. Two exposed cell lines cannot establish population superiority.
[All original arms](../dev/reference_weighted_pseudobulk_20261109/prediction_v3/scores.tsv)
and precedingprediction_v1/v2 outputs remain unchanged.

## Frozen181-line regression replay

All historical single-line queries across both rotations are included. No
post-replay parameter or reference tuning is performed. All181 cases have finite
scores for all three methods, identical within-case genes/windows, and explicit
line identities. Shared reference populations induce dependence; no independent
biological confidence interval or superiority p-value is asserted.

| Pearson summary |Uniform|Uniform +sm5|Uniform +MDL|
|---|---:|---:|---:|
| Mean |.595940|.667955|.668362|
| Median |.614452|.696434|.694386|

Versus unsmoothed uniform, MDL wins178/181 and loses3; mean paired change+.072422,
median+.066858. This establishes mainly the value of smoothing in this replay.
Against the strong sm5 baseline it wins75 and loses106, with mean+.000406 and
median-.006229. The worst change is-.209333(A204_SOFT_TISSUE), best+.236112
(KYSE520_OESOPHAGUS). Harm cannot disappear inside the near-zero mean. The
two locked cases' six relevant scores reproduce exactly in the broader run.

Two independent six-position exhaustive references each enumerate32 partitions
and agree with the dynamic-programming optimum. The documented operation/scale
rounding allowances are7.05e-11 and6.24e-11; objective differences are0. This
checks the partition computation, not Gaussian assumptions or DNA validity.
[Full paired results](../dev/reference_weighted_pseudobulk_20261109/summary_v1/all_paired_profiles.tsv),
[summary and independent checks](../dev/reference_weighted_pseudobulk_20261109/summary_v1/report.json)
and [summary source](../benchmarks/hardening/reference_profile_summary.py)
retain every line. Former evaluations remain regression-only.

## Execution and decision

All five CPU jobs complete0:0:23132423(total weighting),23132548(within-bin),
23132786(MDL locked cases),23133177(frozen181-line replay),23133252(summary).
Inference/prediction jobs use innovation4CPU/16GB/2h, except replay32GB;
summary nslab1CPU/4GB/1h. Actual process peak2348288KiB, replay1m00s.
No GPU, arrays, interactive QOS, installation or frozen-result writes.
Source_v1/v3/v4/v5/v6 respectively execute those five jobs; source_v2 is an
unused intermediate before binding the within-bin design as a manifest input.
Snapshots, neutral-directory imports, wheel/native identities and execution
arguments are retained in the
[job journal](../dev/reference_weighted_pseudobulk_20261109/jobs.json),
[accounting](../dev/reference_weighted_pseudobulk_20261109/accounting.tsv) and
[actual process resources](../dev/reference_weighted_pseudobulk_20261109/resource_summary.json).

Retain implementations as research candidates. No default, release readiness,
scientific threshold, comparator family or historical result changes. The broad
strongest-baseline objective remains unmet; further complexity is not justified
by the two favorable examples or the pooled mean. Biological validation and
uncertainty coverage have not been established by these prediction experiments.
