# RNA-selected smoothing improves the frozen pseudobulk regression

A new RNA-only smoothing selector improves the fixed five-bin baseline across
the frozen181-line Kinker replay:95 wins,33 losses,53 ties, mean paired Pearson
gain.012263 and median.001094. Mean correlation is.680218 versus.667955. The
gain also survives fixed-width and global-RNA-width comparisons. Retain this
as a performance candidate for native-workflow integration, not as an installed
package default, independent validation, DNA-call accuracy or clone-recovery claim.

## Estimand and observation assumptions

The estimand is a regional relative RNA-abundance contrast to the declared
reference. A gene's prepared log ratio is treated as a regional component plus
gene-specific deviation. Neighbouring genes must share enough of that component
for masked prediction to select a useful scale. Absolute error and medians are
robust working choices, not an NB observation law. The Kinker integers were
reconstructed from CPM; no raw-count calibration is asserted.

Regional transcriptional programmes, cis regulation and reference mismatch can
also predict neighbouring genes. RNA prediction alone cannot distinguish those
from a DNA alteration. Independent corresponding DNA at declared resolution is
needed for genomic interpretation. Here the same uncertain relative DepMap DNA
is used only for the already exposed regression assessment. Different cultures,
preselected features and shared references remain limitations. Biological clone
validation and new interval-coverage experiments: not run.

## Fixed algorithm and leakage boundary

[Implementation](../benchmarks/hardening/rna_predictive_smoothing.py) uses the
previously frozen gene log ratios, identities, bins and reference population.
Within each existing bin, genes in their recorded order receive rank-modulo10
folds. Each fold excludes its genes, forms remaining-gene bin medians, and
predicts the excluded genes with within-arm/block running medians. Every eligible
gene is excluded once. Singleton bins contribute no prediction targets, but
remain in the final profile with that support limitation explicit.

The width family is1 and2^k+1(k>=1), through the largest block length. It spans
scales geometrically and includes the fixed five-bin comparator. Choose one
specimen-wide width minimizing summed absolute masked-gene error, then refit
all genes at that width. Numerically tied risk sums choose the smaller width
using an operation/magnitude summation allowance. No DNA-selected grid, width,
reference, uncertainty threshold or per-case DNA winner is used. Every RNA
decision and prediction is saved and hashed before DNA values are parsed.

This is **conditional masked prediction**, not independent confirmation. The
feature universe, normalization, binning and reference had already been estimated
from all source data. The perturbation test checks exclusion at the prepared-RNA
selector boundary only; it does not prove raw preprocessing is leakage-free or
repair the full clone discovery/confirmation problem.

The [initial design](../dev/rna_predictive_smoothing_20261109/DESIGN.md),
[frozen regression continuation](../dev/rna_predictive_smoothing_20261109/REGRESSION_DESIGN.md)
and [post-hoc characterization design](../dev/rna_predictive_smoothing_20261109/CHARACTERIZATION_DESIGN.md)
retain those boundaries. Former evaluations stay regression-only; the candidate
is unchanged after the181-line result.

## Locked cases and complete regression

| Case |Fixed five bins|RNA-selected|RNA-selected width|
|---|---:|---:|---:|
|2313287_STOMACH|.748763|.787832|9|
|42MGBA_CENTRAL_NERVOUS_SYSTEM|.776089|.800024|17|

Both improve, so the prespecified broader replay proceeds. All181 cases and all
baseline arms remain. Predictions use identical DNA-callable genes and1Mb
windows. Widths selected are3(27 lines),5(53),9(56),17(37),33(8). Ties with fixed5
are expected when the procedure selects5; they are not removed from effects.

Five compute fixtures pass: known noiseless width1 optimum; constant-signal
smallest-width tie rule; independent Python-list medians; arm boundary isolation;
and excluded-value perturbations leaving training predictions identical. No
existing test or scientific criterion is changed.

| Method |Mean181-line Pearson|
|---|---:|
|Unsmoothed uniform bin medians|.595940|
|Fixed3|.655166|
|Fixed5|.667955|
|Fixed9|.667943|
|Fixed17|.648592|
|Fixed33|.613130|
|One global width selected from RNA(9)|.667943|
|Global RNA width learned outside each reference group|.667547|
|Per-specimen RNA-selected width|.680218|

Thus the observed mean gain is not explained by substituting one universally
wider fixed window. Against fixed9 the candidate has88 wins,37 losses,56 ties;
mean paired gain.012275, median0. Against the leave-reference-group-out RNA
baseline it has93 wins,38 losses,50 ties; mean+.012672, median+.001518. These
secondary comparisons are descriptive; no new p-value or comparator chosen from
DNA is used to retune the candidate. All362 adaptive/fixed5 scores reproduce
exactly in the characterization. The prior frozen MDL comparator remains visible
in the main replay at mean.668362.

## Shared-reference uncertainty and its limits

Cases are connected if they share any reference or query biological line ID,
including across rotations. This yields8 dependence groups with18,21,24,24,24,
24,24,22 lines. Windows and cells are not replication units. Bootstrap complete
groups, retaining all their lines and using the ratio of summed paired effects
to sampled line count. With20000 draws(seed2026110911), the conditional percentile
95% interval for the line-weighted mean adaptive-minus-fixed5 effect is
[.008356,.016374].

Binomial order-statistic Monte Carlo brackets for the interval endpoints are
[.008263,.008404] and[.016332,.016488], with simultaneous endpoint alpha.0125.
These describe simulation precision of this bootstrap, not biological uncertainty.
The bootstrap assumes independent reference groups and conditions on the selected
feature universe/preprocessing/reference design. Eight groups and exposed DNA
cannot establish prospective coverage, universal superiority or a genomic error
guarantee. Thirty-three harmed lines remain in the full score table.

[Full regression report](../dev/rna_predictive_smoothing_20261109/regression_v1/report.json),
[all case/method scores](../dev/rna_predictive_smoothing_20261109/regression_v1/scores.tsv),
[width ablations, groups and uncertainty](../dev/rna_predictive_smoothing_20261109/characterization_v1/report.json)
and [characterization source](../benchmarks/hardening/rna_smoothing_characterization.py)
retain all evidence. No loss-recall, precision or achieved DNA-call FPR improvement
has been established by Pearson correlation alone.

## Execution and next gate

Three CPU jobs complete0:0:23133691(locked),23134195(all181),23134691(width/group
characterization). Initial/replay innovation4CPU/16GB/2h; characterization
nslab4CPU/8GB/2h. Immutable source_v1/v2/v3, validated noneditable wheel, neutral
working directories, no arrays/shared installs/frozen writes. Process peak649000KiB.
[Commands](../dev/rna_predictive_smoothing_20261109/jobs.json),
[accounting](../dev/rna_predictive_smoothing_20261109/accounting.tsv) and
[actual resources](../dev/rna_predictive_smoothing_20261109/resource_summary.json).

The installed package and its native exports have not yet adopted this selector.
The next gate is native-workflow integration and direct comparison at the same
measurement/population, followed by task-specific detection/error assessment and
appropriate independent evaluation. Do not transfer this correlation gain into
clone confidence, spatial genotype assignments or count-model calibration flags.
