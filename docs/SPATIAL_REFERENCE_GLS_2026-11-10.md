# Fit mixed spatial spots using reference-estimated variance

## Result

Job 23185898 completed in 28 seconds on CPU node003. All 128 simulated
specimens, 1,024 method rows and 112 planned comparisons completed.
The independent rational covariance, simplex solutions and coordinate
permutation checks passed. The process peaked at 521,760 KiB. Slurm
reported zero batch RSS.

| Reads | Count and state law | Reference GLS MAE | Averaged multinomial mixture MAE | Averaged least-squares MAE |
|---:|---|---:|---:|---:|
| 750 | Multinomial | 0.122383 | 0.121686 | 0.151706 |
| 750 | Overdispersion only | 0.155048 | 0.162779 | 0.209381 |
| 750 | State only | 0.253586 | 0.253408 | 0.250987 |
| 750 | State and overdispersion | 0.261309 | 0.259209 | 0.262629 |
| 3000 | Multinomial | 0.074331 | 0.074234 | 0.094115 |
| 3000 | Overdispersion only | 0.119813 | 0.148083 | 0.212252 |
| 3000 | State only | 0.259522 | 0.258338 | 0.253742 |
| 3000 | State and overdispersion | 0.267495 | 0.261951 | 0.256719 |

At 3,000 reads with overdispersion alone, both error endpoints improve
in all 16 specimens against both strong averaged comparators. Regional
profile L1 is 0.014084 versus 0.017685 for the multinomial mixture and
0.026253 for least squares. At 750 reads with overdispersion alone, the
point gains against the averaged multinomial mixture fail the adjusted
sign criterion. Mean reference estimates of dispersion are 0.10116 and
0.09748 at the two depths under the Gamma-noise law, whose input is 0.1.
Under ordinary multinomial sampling the mean estimates are 0.00060 and
0.00015.

The fit passes 89 of 112 comparisons and fails the full gate. With the
transcriptional field present, neither depth improves regional profile
error against averaged least squares, even when the independent reference
captures overdispersion. This factorial result shows a specific benefit
of variance weighting for extra count noise and an unresolved shift in
expected RNA under the simulated state law. Separating the genomic
interpretation from that shift would require an independent state
measurement or matched DNA.

The estimator remains experimental. The simulation uses known pure
profiles and zero biological specimens. All frozen outputs, scientific
criteria, and earlier negative gates were preserved.

[All 112 contrasts](../dev/spatial_reference_gls_20261110/study_v1/report.json),
[all 1,024 method rows](../dev/spatial_reference_gls_20261110/study_v1/scores.tsv),
[mean errors by condition](../dev/spatial_reference_gls_20261110/study_v1/means.tsv),
[numerical checks](../dev/spatial_reference_gls_20261110/study_v1/numerics.json),
[completed manifest](../dev/spatial_reference_gls_20261110/study_v1/run_manifest.json),
[job accounting](../dev/spatial_reference_gls_20261110/accounting_complete.tsv).

## Design and execution

Job 23185898 was submitted on nslab CPU with 1 CPU, 8 GB RAM, and a
one-hour limit. Driver and launcher syntax checks passed. The immutable
snapshot preserves all 959 parent files.

The preceding spatial experiment combined transcriptional state and
overdispersed counts. This new panel crosses those two effects. All
methods see the same fixed profiles, spots, counts and spatial graph
within each simulation. Across conditions, the sixteen seed indices
share the same underlying composition within each depth. There are
128 scenario-specific simulations and zero biological specimens.

The new estimator uses 100 independent reference-only spots per
seed, depth and noise setting to estimate variance above ordinary
multinomial sampling. It whitens the compositions using that working
covariance, fits the same simplex least-squares model, and applies
the same spatial averaging. The fixed comparison covers both
contribution and profile error against all seven prior methods.
The count-variance correction targets noise. A CNA-shaped state can
still shift the expected RNA composition.

[Frozen design](../dev/spatial_reference_gls_20261110/DESIGN.md),
[job command](../dev/spatial_reference_gls_20261110/jobs.json),
[source snapshot](../dev/spatial_reference_gls_20261110/snapshot.json),
[driver](../benchmarks/hardening/spatial_reference_gls.py).
