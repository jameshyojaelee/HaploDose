# Estimate RNA contributions in mixed spatial spots

## Fixed-profile result

Job 23185815 completed in 42 seconds on CPU node003. All 64 specimens
and 40 paired contrasts are present. The analytic, permutation,
identifiability and observational-equivalence fixtures passed. The
process peaked at 520,692 KiB and Slurm recorded 258,868 KiB.

| Reads per spot | Observation law | Mixture contribution MAE | Native spatial NB MAE | Neighbour least-squares MAE |
|---:|---|---:|---:|---:|
| 750 | Matched | 0.198816 | 0.215984 | 0.147479 |
| 750 | State and overdispersion | 0.322715 | 0.340600 | 0.250072 |
| 3000 | Matched | 0.120867 | 0.218462 | 0.098357 |
| 3000 | State and overdispersion | 0.313903 | 0.343773 | 0.260924 |

The mixture improves both endpoints against the native spatial NB
assignments in all four point comparisons. Only the matched 3,000-read
condition passes both adjusted sign tests against that comparator.
Neighbour-averaged least squares leads the unsmoothed mixture on both
mean errors in all four conditions. Reject the unsmoothed method under
the declared joint gate. The synthetic stress represents a competing
RNA explanation. It provides zero biological specimens.

[All 40 comparisons](../dev/spatial_fixed_profile_mixture_20261110/study_v1/report.json),
[all 384 score rows](../dev/spatial_fixed_profile_mixture_20261110/study_v1/scores.tsv),
[numerical fixtures](../dev/spatial_fixed_profile_mixture_20261110/study_v1/numerics.json),
[complete manifest](../dev/spatial_fixed_profile_mixture_20261110/study_v1/run_manifest.json).

## Matched spatial averaging

Job 23185834 completed in 15 seconds on nslab CPU node003 with 1 CPU,
8 GB RAM, and a 30-minute limit. Process peak RSS was 446,512 KiB.
Slurm reported zero RSS for the batch step. All 64 archived simulations
match their hashes. All 384 prior score rows and 1,152 error endpoints
reproduce exactly. The graph reproduces each stored neighbour
least-squares prediction exactly.

| Reads per spot | Observation law | Averaged mixture MAE | Averaged least-squares MAE | Averaged mixture profile L1 | Averaged least-squares profile L1 |
|---:|---|---:|---:|---:|---:|
| 750 | Matched | 0.120967 | 0.147479 | 0.012953 | 0.016176 |
| 750 | State and overdispersion | 0.252808 | 0.250072 | 0.033484 | 0.031857 |
| 3000 | Matched | 0.076529 | 0.098357 | 0.008528 | 0.010550 |
| 3000 | State and overdispersion | 0.254804 | 0.260924 | 0.038738 | 0.037857 |

The averaged mixture passes 44 of 48 comparisons. It improves both
endpoints in every specimen against native spatial NB assignment,
including both stressed conditions. It also passes both endpoints at
both depths against averaged least squares under the matched law. The
four failed contrasts are the two endpoints at both depths against
averaged least squares under stress. At 3,000 reads its point MAE is
lower by 0.006120, with 9 of 16 certified wins and adjusted p = 1.
Its profile error is higher under stress at both depths. Boundary MAE
is lower than averaged least squares under the matched law and higher
under stress at both depths.

This identifies a benefit of explicit mixtures with matched smoothing
when the pure RNA profiles are known and the likelihood fits the count
law. The combined stress does not isolate transcription from
overdispersion. The full gate fails and the research implementation
is withheld. No biological specimens, DNA correspondence or new clone
calls enter either study. The 5% specimen false-split criterion, power,
and frozen results are preserved.

[Matched design](../dev/spatial_matched_mixture_smoothing_20261110/DESIGN.md),
[matched job command](../dev/spatial_matched_mixture_smoothing_20261110/jobs.json),
[all 48 comparisons](../dev/spatial_matched_mixture_smoothing_20261110/study_v1/report.json),
[all 448 score rows](../dev/spatial_matched_mixture_smoothing_20261110/study_v1/scores.tsv),
[matched manifest](../dev/spatial_matched_mixture_smoothing_20261110/study_v1/run_manifest.json),
[job accounting](../dev/spatial_fixed_profile_mixture_20261110/accounting_complete.tsv).

## Original submission and design

Submitted job 23185815 on the CPU partition with innovation, 1 CPU,
8 GB RAM, and a one-hour limit. The driver and launcher passed syntax
checks. The immutable source preserves all 955 parent files.

The test supplies the same three fixed normalized RNA profiles to every
estimator and uses a mixture of their expected contributions for each
spot. It compares the fitted simplex likelihood against direct and
neighbour-averaged least squares, unspatial categorical responsibilities,
and both multinomial and NB inputs to native spatial smoothing. All six
methods use identical spots and fixed profiles. The outcomes are RNA
contribution error and regional RNA composition error. Converting to
cell proportions would require transcript output per cell.

All 64 simulated specimens are prespecified across two read depths and
two observation laws. The stressed law includes a spatially structured
transcriptional effect over the loss region and Gamma rate noise. The
exact observational-equivalence fixture assigns the same probability
vector to a transcriptional state and a genomic mixture. Distinguishing
those explanations requires independent measurements. Numerical checks
use analytic interior, boundary and pure-spot solutions plus a gap and
curvature-based error allowance. The joint gate uses all 40 paired
contrasts. The original clone false-split and power criteria still apply.

[Frozen design](../dev/spatial_fixed_profile_mixture_20261110/DESIGN.md),
[job command](../dev/spatial_fixed_profile_mixture_20261110/jobs.json),
[snapshot](../dev/spatial_fixed_profile_mixture_20261110/snapshot.json),
[driver](../benchmarks/hardening/spatial_fixed_profile_mixture.py).
