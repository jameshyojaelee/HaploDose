# Reliability weighting for registered spatial reference RNA

Reference reliability weighting reduces both errors relative to the
unshrunk spatial correction in all four state-absent conditions. It
passes 113/128 contrasts, up from 110/128, and retains all 64
state-present passes. The unchanged joint gate fails. Against the
variance-only GLS control, none of the eight state-absent endpoint
contrasts passes. Under shared state with overdispersion, shrinkage
also increases both mean errors relative to the unshrunk correction.
The candidate remains a research result.

## Mechanism and decision

The two reference groups use the same registered band positions and
disjoint independently simulated spots at equal sequencing depth. Their
cross product estimates shared field strength in the fixed inverse
baseline-proportion metric. The squared norm of their pooled mean
includes reference estimation noise. Their ratio, clipped to [0, 1],
shrinks the registered correction toward the fixed baseline. A single
coefficient applies to all bands and genomic blocks. This regularizes
the reference measurement and leaves the query fitting, neighbourhood
averaging and profile definitions intact. The ratio itself has no
calibration guarantee. The inference assumes the state field transfers
from the independently registered reference to the query.

| Reads | Observation law | Median reference weight | Mean gain versus unshrunk in contribution MAE | Mean gain versus unshrunk in profile L1 | Passes under original gate |
|---:|---|---:|---:|---:|---:|
| 750 | Multinomial | 0.173 | 0.003956 | 0.000575 | 11/16 |
| 750 | Overdispersion | 0.081 | 0.006175 | 0.000873 | 12/16 |
| 750 | State | 0.943 | 0.000869 | −0.000048 | 16/16 |
| 750 | State and overdispersion | 0.840 | −0.004498 | −0.000838 | 16/16 |
| 3,000 | Multinomial | 0.006 | 0.004502 | 0.000493 | 12/16 |
| 3,000 | Overdispersion | 0.000 | 0.012639 | 0.001551 | 14/16 |
| 3,000 | State | 0.983 | −0.000242 | −0.000083 | 16/16 |
| 3,000 | State and overdispersion | 0.892 | −0.003775 | −0.000890 | 16/16 |

The first two numeric gain columns are descriptive paired differences
against the unshrunk correction. The 128 primary contrasts retain the
original eight comparators, two endpoints, both depths and four laws,
with sign tests, numerical allowances and the original Bonferroni
criterion. The variance-only comparator has 3–8 certified wins for the
candidate out of 16 in state-absent conditions. The largest improvement
in mean error over that comparator is 0.000800 for contribution MAE,
at 750 reads under overdispersion. These exposed simulations give
mechanistic evidence about reference noise, with zero biological
specimens and no change in the spatial inference release.

## Verification and resources

All 128 pooled multipliers reproduce their original values exactly
before shrinkage. A rational test gives weight 3/4, and swapping the
halves preserves it. Opposed estimates return zero and identical
nonzero estimates return one. All 128 candidate reference files are
hash-locked before query counts and truth are loaded. The complete
score table preserves all 1,280 old rows exactly. Job 23188033 completes
0:0 in 29 seconds on nslab CPU node003 with 1 CPU, 8 GB and a one-hour
limit. The process records 452,420 KiB peak RSS. The job's short duration
leaves Slurm batch RSS at zero.

[Design](../dev/spatial_reference_shrinkage_20261110/DESIGN.md),
[all contrasts](../dev/spatial_reference_shrinkage_20261110/run_23188033/results/report.json),
[score table](../dev/spatial_reference_shrinkage_20261110/run_23188033/results/scores.tsv),
[locked references](../dev/spatial_reference_shrinkage_20261110/run_23188033/results/REFERENCES_LOCKED.json),
[manifest](../dev/spatial_reference_shrinkage_20261110/run_23188033/results/run_manifest.json),
[job command](../dev/spatial_reference_shrinkage_20261110/jobs.json),
[resource log](../dev/spatial_reference_shrinkage_20261110/run_23188033/resources.txt).
