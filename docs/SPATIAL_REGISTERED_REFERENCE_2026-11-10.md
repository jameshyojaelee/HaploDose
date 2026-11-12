# Estimate spatial transcriptional state from registered reference RNA

## Result

Job 23185965 completed in 20 seconds on CPU node003. All 128 independent
reference calibrations were locked before query counts or truth were
decoded. All 1,024 preceding score rows and 3,072 error endpoints
reproduced exactly. Rational multiplier and component transformation
fixtures passed. The process peaked at 447,068 KiB. Slurm reported
zero batch RSS.

| Reads | Observation law | Variance-only contribution MAE | Registered reference MAE | Reversed reference MAE | Variance-only profile L1 | Registered reference L1 |
|---:|---|---:|---:|---:|---:|---:|
| 750 | State only | 0.253586 | 0.126504 | 0.339626 | 0.033789 | 0.013227 |
| 750 | State and overdispersion | 0.261309 | 0.157659 | 0.348130 | 0.034282 | 0.016958 |
| 3000 | State only | 0.259522 | 0.077728 | 0.349261 | 0.039382 | 0.009086 |
| 3000 | State and overdispersion | 0.267495 | 0.134521 | 0.360413 | 0.039545 | 0.015617 |

The correctly registered reference improves both endpoints in all 16
replicates against every comparator in every state-present stratum.
All 64 contrasts in those strata pass the family-adjusted sign tests.
Reversing the orientation worsens both mean endpoints relative to
the variance-only fit in all four state-present strata. These results
depend on the simulated shared reference-to-tumour state response.

The full gate passes 110 of 128 contrasts. All 18 failures occur when
the state is absent. In those strata, the registered correction has
higher mean contribution and profile error than variance-only GLS at
both depths. The experiment withholds promotion. It distinguishes the
benefit of independently measured matching state from the cost of
injecting reference noise or registering a different state field.
All prior gates, original clone criteria and frozen results remain
preserved. Biological sample size is zero.

[All 128 contrasts](../dev/spatial_registered_reference_20261110/study_v1/report.json),
[all 1,280 score rows](../dev/spatial_registered_reference_20261110/study_v1/scores.tsv),
[reference lock](../dev/spatial_registered_reference_20261110/study_v1/REFERENCES_LOCKED.json),
[completed manifest](../dev/spatial_registered_reference_20261110/study_v1/run_manifest.json),
[job accounting](../dev/spatial_registered_reference_20261110/accounting_complete.tsv).

## Candidate independent spatial evidence

The authors' archive lists the human colorectal Figure 4 DNA section
`human_colon_cancer_4_dna_200114_13` and RNA section
`human_colon_cancer_4_rna_200102_06`. The Figure 3 DNA and scWGS entries
come from a different sampled region of the same tumour. The Figure 4
paper describes co-registered serial sections. This provides a lead
for independent regional DNA compatibility, with one human tumour,
registration uncertainty, and no one-to-one DNA/RNA spot correspondence.
Availability of a separately registered reference population with a
transferrable transcriptional response is unverified. Measurements
remain unopened and the cohort is not yet admitted for evaluation.
[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9301586/),
[official sample manifest](https://singlecell.broadinstitute.org/single_cell/study/SCP1278/spatial-genomics-enables-multi-modal-study-of-clonal-heterogeneity-in-tissues),
[authors' archived code](https://doi.org/10.5281/zenodo.5553305).

## Original design and execution

Job 23185965 was submitted on innovation CPU node003 with 1 CPU, 8 GB
RAM, and a one-hour limit. Driver and launcher syntax checks passed.
The immutable snapshot preserves all 961 parent files.

The preceding variance-weighted fit still lost profile accuracy when
a spatial transcriptional field overlapped the simulated loss. This
experiment supplies independently generated pure reference spots on
the same coordinate grid. It estimates regional changes solely from
reference counts and locks all estimates before reading query counts.
The state-corrected simplex fit adjusts each fixed component profile
and recovers contributions on the original RNA scale before averaging.

The identifying assumption is a shared multiplicative transcriptional
response at each registered position. An additional control reverses
the reference orientation while retaining the same reference reads.
All 128 prior queries, eight comparators, exact score controls and the
fixed joint decision remain. This simulation has zero biological
specimens. Transfer in tissue would require independent registration
and evidence that the reference captures the relevant state.

[Frozen design](../dev/spatial_registered_reference_20261110/DESIGN.md),
[job command](../dev/spatial_registered_reference_20261110/jobs.json),
[source snapshot](../dev/spatial_registered_reference_20261110/snapshot.json),
[driver](../benchmarks/hardening/spatial_registered_reference.py).
