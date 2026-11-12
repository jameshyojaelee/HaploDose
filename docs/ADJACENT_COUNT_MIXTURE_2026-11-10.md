# Adjacent-region count mixtures for individual-cell RNA profiles

At 3,000 retained molecules and without extra RNA state, the fitted
mixture retains 94.5% of the generated between-population contrast.
Regional shrinkage retains 62.9% on the same fresh draws. The mixture
loses its minority-region advantage under continuous RNA state and has
little whole-profile benefit over pseudobulk at low depth. It passes
24 of 48 prespecified comparisons. The joint gate fails.

## Measurement and test

An individual cell's relative regional RNA log odds is the target. The
candidate learns shared compositions from two other cell folds. It selects
one of two adjacent partitions and two or three mixture components by
held-out cell count likelihood. When neither partition improves likelihood
over one component, the output is leave-one-cell-out pseudobulk. The
active mixture estimates a cell's composition from its own counts and
the two training folds. It makes no genomic membership decision.

The fresh simulation uses 48 independent seeds per condition, with depths
of 300 and 3,000 paired within seed. Minority prevalence is 20%. Event
starts are randomized within chromosome. The stress condition introduces
continuous RNA state across two chromosomes, including the minority
gain region. The scorer compares whole-profile and minority-region RNA
mean squared errors. Six comparators include the same mixture procedure
without grouping and a phase-averaged adjacent count ratio. The declared
family uses simultaneous one-sided exact binomial bounds with Bonferroni.
The lower bound on win probability must exceed 0.5 for every contrast.

| Condition | Depth | Whole error, mixture / pseudobulk | Minority error, mixture / individual | Generated contrast retained, mixture / regional shrinkage |
|---|---:|---:|---:|---:|
| Two CNA populations | 300 | 0.004157 / 0.004192 | 0.075143 / 0.115943 | 10.5% / 16.2% |
| Two CNA populations | 3,000 | 0.000357 / 0.003470 | 0.002576 / 0.011098 | 94.5% / 62.9% |
| Two CNA populations and state | 300 | 0.010753 / 0.016213 | 0.085613 / 0.113415 | 12.2% / 25.0% |
| Two CNA populations and state | 3,000 | 0.006419 / 0.015500 | 0.069538 / 0.010748 | 21.8% / 73.3% |

In the low-depth two-population condition, 112 of 144 fits fell back to
pseudobulk. At the specimen level, 26 of 48 whole-profile results tied
pseudobulk exactly, and ten improved. At high depth without added state,
the mixture wins against the individual ratio in 47/48 minority-region
comparisons and against pseudobulk in 48/48. At high depth with state,
the mixture wins 0/48 on minority-region error against the individual
ratio. It selects three components in 142/144 fits. A model optimized
for the full RNA composition can assign mixture capacity to continuous
transcriptional variation. This experiment does not identify which
fitted component corresponds to a genomic group. The ungrouped mixture
also has large minority error under this stress, and no joint contrast
against it passes. In homogeneous controls the grouped mixture falls back
in 143/144 and 141/144 fits at the two depths. These are prediction choices,
outside the specimen false-split endpoint.

## Validation and execution

The rational-posterior, chromosome-boundary, count-conservation and
held-out-fold checks pass. All 1,152 previous score rows reproduce exactly.
Doubling the iteration budget in 18 prespecified strata-by-seed fixtures
changes no selected mixture structure. The largest adjacent-mixture error
changes are 0.00000174 whole profile and 0.00000915 minority. This assesses
iteration sensitivity on those fixtures. The full fitted mixture remains
an experimental predictor.

Job 23188288 completed `0:0` in 27 seconds on `ne1dc6-010` under `nslab` CPU
with 1 CPU, 8 GB and 1 hour. The study used 82,176 KiB peak process RSS and
wrote 2,016 score rows. There are zero biological specimens. Production
source, frozen results, and all existing scientific criteria were preserved.

[Design](../dev/adjacent_count_mixture_20261110/DESIGN.md),
[source](../dev/adjacent_count_mixture_20261110/study.py),
[submission](../dev/adjacent_count_mixture_20261110/submission.json),
[accounting](../dev/adjacent_count_mixture_20261110/accounting.tsv),
[summary](../dev/adjacent_count_mixture_20261110/run_23188288/results/report.json),
[comparisons](../dev/adjacent_count_mixture_20261110/run_23188288/results/comparisons.tsv),
[scores](../dev/adjacent_count_mixture_20261110/run_23188288/results/scores.tsv),
[fit choices](../dev/adjacent_count_mixture_20261110/run_23188288/results/selected_models.json),
[iteration checks](../dev/adjacent_count_mixture_20261110/run_23188288/results/iterations.tsv).
