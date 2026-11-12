# Reference weighting, gene aggregation and membership information

The scientific goal remains supported RNA-derived CNA profiles and clone
assignments, with uncertainty and independent evidence at their claimed
resolution. This round rejects two proposed bulk repairs and separates intrinsic
membership overlap from algorithmic loss. It changes no inference default.

## Bulk mechanism test

The experiment crosses pooled-depth versus equal-donor reference rates with
summed-count versus median-gene-log-ratio aggregation. It reconstructs the same
nine raw-STAR development fits: 20 normal donors, 24 calibration patients and
24 previously inspected reserved patients, with eight reference resampling fits.
No cohort is newly designated independent evaluation. Biological-unit paired
effects use patients; shared reference variation is represented by the frozen
donor refits, not by treating bins as independent replicates.

The estimand is the paired change in DNA-profile agreement or loss recall on
common callable windows. All profiles use the same RNA-only partition operator,
penalty, noise scale, arm/gap breaks and minimum three windows, with their own
boundaries. Existing calibration patients supply loss thresholds. Reserved
outcomes do not select thresholds, references, genes or boundaries.

Equal-donor arithmetic gene rates isolate depth weighting, not the naive
baseline's mean-log(CPM+1) construction. Median aggregation changes both the
elementary observation unit and allocation of the half-count offset. It does
not isolate an offset-only mechanism and is not assigned count-model uncertainty.
These distinctions prevent attributing any difference to the wrong intervention.

All ten earlier comparators and 36 endpoint comparisons remain. Three new
methods and six paired contrasts expand the family to 54 endpoints, with
20,000 nested resamples, Holm tests and simultaneous Bonferroni basic bootstrap
intervals. Eight reference draws and finite tail resolution limit interval
interpretation; these approximate intervals have not established coverage.

| Contrast | Paired Pearson effect | Simultaneous interval |
|---|---:|---:|
| Equal versus pooled reference, summed counts | −0.00036 | [−0.00806, 0.00536] |
| Equal versus pooled reference, median genes | −0.00207 | [−0.00924, 0.00610] |
| Median versus summed aggregation, pooled reference | 0.01464 | [−0.03507, 0.05227] |
| Median versus summed aggregation, equal reference | 0.01293 | [−0.03677, 0.05013] |
| Pooled-reference median versus matched naive | −0.04140 | [−0.06276, −0.02147] |
| Equal-reference median versus matched naive | −0.04347 | [−0.06530, −0.02168] |

Correlations are evaluable in 23/24 patients; loss recall in 13/24. The two
median-versus-naive Pearson contrasts have Holm p=0.00270. No new contrast
shows a simultaneous benefit in Spearman or loss recall. Median variants have
mean recall 0.3140, with loss FPR 0.0464 and 0.0480; naive partition has recall
0.3997 and FPR 0.0412. All retain full classified-truth coverage on common
support; undefined precision is not filled with zero or silently discarded.

**Decision:** neither reference reweighting nor median aggregation is supported
as a repair here. Absence of a significant effect is not equivalence. The
remaining naive-reference transformation differs, and specimen correspondence,
DNA truth uncertainty and unmodeled expression programmes remain possible limits.
The current evidence does not justify adding another default or claiming a
special advantage of the more complex estimator.

## Membership information diagnostic

The completed shared-membership experiment has the same state-calibrated
false-split and detection counts for control and candidate: 12/200 and 88/200.
Its simultaneous upper false-split bound is 12.93%, above the unchanged 5%
criterion. Only 1/200 alternatives has exactly two raw tumour groups. The
candidate is not promoted despite improving the five-case recovery panel to 4/5.

For a secondary diagnostic, the code was frozen before reading the completed
null summary. It reads saved alternatives after fitting, without changing
inference. Under the generator's overlapping Gaussian state amplitudes, equal
independent class priors and a known CNA effect, observing noiseless amplitude
gives expected accuracy Phi(delta/(2 sigma))=0.686714. Independent integration
agrees; this is conditional on that observation model and exchangeable cell
identities. It is not a guarantee for individual specimens, a fixed-balanced
sample theorem, or proof that a difference is genomic rather than transcriptional.

The 200 saved alternatives give mean oracle accuracy 0.684892. A separate
fixed-balanced-count ranking oracle gives 0.685033. One-to-one matching of raw
groups to truth gives 0.438292 for control and 0.439617 for the exact candidate;
extra groups and unassigned cells receive no credit. Independent fixtures check
perfect, oversplit, half-assigned and wholly unassigned cases. Candidate coverage
is 99.61%; imposing the specimen gate as abstention gives 43.81% mean coverage
and 19.10% correct assignments across all cells. Reported specimen percentiles
are variability summaries, not confidence intervals.

[Inference] Both overlap and algorithmic subdivision contribute. More precise
confirmation alone did not close the gap. Tracing initialization and merging is
the next bounded mechanism test; oracle labels remain excluded from inference.
Same-cell DNA/RNA and a sufficiently replicated biological null panel are still
missing for the intended biological guarantees.

## Execution and preservation

Job22963075 completed the factorial in 21m32s; process peak RSS was 680 MiB.
Audit22963110 completed in 15s and verifies 109 input hashes, 382 baseline files,
65 unchanged package/test/native files, ten unchanged earlier comparator outputs
and all 54 paired effects. Independent scalar calculations, gene permutation,
equal-reference invariance to donor depth, count conservation and exact pooled
control reconstruction pass. The numerical allowance derives from float64
accumulation and the independent reference, not a scientific decision cutoff.

Job22963381 completed the membership diagnostic in 27s (137 MiB process peak
RSS). Joint-null audit22963491 additionally reconciles all 45 endpoint counts,
nine paired discordances, 3,600 decisions and independent exact confidence bounds.
The original audit and all failed or negative experiments remain preserved.

New implementation: `benchmarks/hardening/aggregation_reference.py`,
`aggregation_reference_panel.py`, `aggregation_reference_audit.py` and
`membership_information_reference.py`. Immutable sources, designs, parameters,
hashes, wheel identity and scheduler records are under
`dev/aggregation_reference_20261108`. No production dependency or package source
changed. The unchanged full suite was not rerun; last evidence remains
319 Python passes, one GPU-only skip and four Rust passes in job22855685.
