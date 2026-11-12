# Cross-fitted regional shrinkage for individual-cell RNA profiles

The candidate passes 20 of 24 prespecified comparisons. It reduces
minority-region RNA error versus leave-one-cell-out pseudobulk in all four
heterogeneous strata, with 48 wins in 48 simulations per stratum. It loses
to pseudobulk on the whole profile at low depth without extra RNA state. At
higher depth it loses to the individual-cell ratio on minority regions in
both heterogeneous scenarios. The joint gate fails and the estimator stays
outside the production package.

## Design and what was measured

The November 10 review documented a failed individual-cell transfer test
against exposed line-level DNA. Here the target is an individual cell's
regional relative RNA log odds under a known conditional count law. A beta
working prior shrinks each regional probability towards a mean estimated
from another cell fold. Its strength uses the difference between within-cell
factorial moments and cross-cell count products to separate technical
binomial variance from between-cell variation. The latter includes dosage
and transcriptional variation. An individual-cell ratio, a common shrinkage
strength, and leave-one-cell-out specimen pseudobulk provide matched controls.

Each regime has 48 independently generated specimens. The depths of 300 and
3,000 retained molecules are paired across replicate seeds. Each query has
120 cells and 24 regions with matched support and an independent normal count
pool. The heterogeneous conditions add a minority gain and loss to 20% of
cells. Continuous state overlaps the minority gain in the stress condition.
The primary endpoints are whole-profile and minority-region mean squared
error. The frozen 24-comparison family uses exact one-sided binomial lower
bounds with Bonferroni adjustment. Each lower bound must exceed 0.5. Cell and
region counts are nested within a simulation specimen.

## Results

Values are mean squared relative RNA log-odds errors across simulations.

| Condition | Depth | Regional whole | Pseudobulk whole | Regional minority | Pseudobulk minority | Individual minority |
|---|---:|---:|---:|---:|---:|---:|
| Homogeneous | 300 | 0.002551 | 0.000916 | | | |
| Homogeneous | 3,000 | 0.000244 | 0.000079 | | | |
| Two CNA populations | 300 | 0.005176 | 0.004133 | 0.067301 | 0.086339 | 0.114190 |
| Two CNA populations | 3,000 | 0.001273 | 0.003485 | 0.014862 | 0.084446 | 0.010577 |
| Two CNA populations and state | 300 | 0.011305 | 0.016196 | 0.068168 | 0.117089 | 0.108649 |
| Two CNA populations and state | 3,000 | 0.002325 | 0.015545 | 0.014456 | 0.115444 | 0.010325 |

At low depth in the two-population condition without extra state, the
candidate wins on the whole profile in 0/48 specimens versus pseudobulk and
3/48 versus global shrinkage. The homogeneous control also favors pseudobulk.
These are expected costs of estimating regional heterogeneity from noisy
counts. At higher depth the candidate wins on the minority endpoint against
the individual ratio in 1/48 specimens without state and 2/48 with state.
The fitted prior pulls minority cells towards the specimen mean. In the
no-state condition it retains a mean 15.5% of the generated group contrast at
low depth and 65.7% at high depth. The individual ratio retains 98.4% and
100.9% respectively, with much larger whole-profile error. This identifies
two practical failure modes for this shrinkage rule.

With continuous state, regional shrinkage improves RNA prediction while
giving higher state-free response error than pseudobulk at both depths.
The observed counts can support the regional RNA effect even when that effect
mixes transcriptional and genomic responses. There are zero biological
specimens and no membership, calibration, or false-split claim.

## Execution and preservation

Slurm job 23188235 completed `0:0` in 23 s on `ne1dc6-010` under `nslab` CPU
with 1 CPU, 8 GB and 1 h. The Python study used 79,872 KiB peak process RSS.
The exact moment/posterior and limiting-case checks pass. A held-out cell
perturbation and a row permutation confirm fold isolation and output identity.
The 1,152 score rows and all 24 contrasts are saved. The job writes only to
its new run directory. Production source and frozen results remain preserved.

[Frozen design](../dev/crossfit_regional_profile_20261110/DESIGN.md),
[source](../dev/crossfit_regional_profile_20261110/study.py),
[submission](../dev/crossfit_regional_profile_20261110/submission.json),
[report](../dev/crossfit_regional_profile_20261110/run_23188235/results/report.json),
[paired contrasts](../dev/crossfit_regional_profile_20261110/run_23188235/results/comparisons.tsv),
[all scores](../dev/crossfit_regional_profile_20261110/run_23188235/results/scores.tsv),
[accounting](../dev/crossfit_regional_profile_20261110/accounting.tsv).
