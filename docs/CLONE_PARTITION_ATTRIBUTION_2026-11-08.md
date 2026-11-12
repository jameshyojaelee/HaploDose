# Saved assignments show fragmentation and membership mixing

All 800 saved state-null and two-clone evaluation cases from the paired fresh
joint experiment were inspected, without refitting or selecting cases. Job
23000251 completes successfully in 30 seconds, with 114428 KiB peak process RSS.
Chunk and assignment hashes are checked. All 800 assignment matchings agree
with independent exhaustive enumeration. The design precedes execution in
`dev/clone_partition_attribution_20261108/DESIGN.md`.

## Descriptive estimands and results

These are raw final assignments, before the specimen-level calibration gate.
Correctness denominators include all 600 simulated cells, retaining unassigned
and reference-like cells. One-to-one matching maps at most one output group to
each true clone. Oracle many-to-one matching lets every eligible output group
choose its majority true label, using unavailable simulation truth. It measures
the maximum correctness obtainable by simply relabelling or coarsening these
fixed groups. It is not a usable assignment method or a successful-recovery test.

| Arm and condition | Mean groups | Coverage | One-to-one correctness | Oracle many-to-one correctness | Within-truth state variance explained |
|---|---:|---:|---:|---:|---:|
| Proposal merge, state null | 1.785 | 93.23% | 68.85% | 93.23% | 39.81% |
| Full merge, state null | 1.790 | 93.08% | 68.81% | 93.08% | 39.68% |
| Proposal merge, two clones | 3.135 | 97.34% | 46.39% | 64.23% | 66.68% |
| Full merge, two clones | 3.150 | 97.34% | 46.21% | 64.24% | 66.76% |

Each row summarizes 200 conditional simulated specimens. No biological specimens
are added. State association is the descriptive variance explained by final
groups after centering the known latent state separately within each true CNA
population. It includes all MAP labels, including reference-like assignments;
adaptive selection and group counts affect this association. It does not identify
a successful nuisance correction. No new significance tests, cutoffs or claims of
biological calibration are introduced.

In two-clone cases the full-merge arm has mean covered-weighted conditional truth
entropy 0.8618 bits. Its oracle coarsening gap is 18.03 percentage points, yet
oracle correctness remains only 64.24%. Thus both fragmentation and within-group
mixing are present. These fixed assignments cannot become accurate merely by
choosing a more favorable merge. The nearly identical two arms also indicate that
preserving full memberships at merge did not resolve the main information loss.

## Decision

Retain the original three thinning failures and failed false-split/power gates.
Do not promote the joint full-membership procedure. Further method work must
target information used to form memberships and independently constrain regional
state variation; a merge-only intervention is insufficient for these saved groups.
The earlier conditional information oracle remains relevant. This diagnosis does
not imply that the public RNA workflow can separate DNA clones from arbitrary
transcriptional state, or that more compute can replace independent DNA evidence.

Reproducible records are in `dev/clone_partition_attribution_20261108/diagnostic_v1`:
`cases.tsv`, `descriptive_summary.tsv`, `report.json` and `run_manifest.json`.
The script is `benchmarks/hardening/clone_partition_attribution.py`. Every original
fit, assignment, calibration decision and scientific criterion is preserved.
