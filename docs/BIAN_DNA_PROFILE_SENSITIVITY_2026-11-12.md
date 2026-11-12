# Bian DNA profiles: fixed measurement-sensitivity analysis

The descriptive DNA normalization and processing-sensitivity analysis is now
implemented and its known-answer fixtures pass on compute. The full analysis
now completes after exact native reconstruction. No RNA/DNA performance has been scored,
no genomic truth has been admitted, and no scientific method has been promoted.
This is a prerequisite already specified by the full-library measurement plan,
not a replacement for the requested genomic accuracy and clone-recovery goal.

## Completed full result and admission decision

Native23233877 completes0:0, all274 before/after identities and all24 exact
count comparisons pass. Full sensitivity **23235024 completes0:0**,
nslab/nslab CPU1/8GiB/2h,ne1dc5-013. Slurm elapsed1m56s; timed1m54.76s,
process peak97948KiB and sampled Slurm peak80908KiB. Existing environment reports
Python3.11.18,NumPy2.4.6,Pandas2.3.3; backend CPU NumPy/Pandas/Decimal.
Dataframe-fragmentation performance warnings remain in the original log.
No inference kernel, GPU, new dependency, batch submission or public patch.

All **312** before/after source identities match. Actual native-input hashes
match before/after. All integer nesting/conservation gates pass. Every one of
the **204** full/split numerical checks passes; maximum Decimal discrepancy
8.881784197001252e-16, declared arithmetic budgets1.9191e-14–3.2015e-14.
All12 settings,17 vectors each,20 processing comparisons and projections onto
1630/1166 frozen RNA regions remain. The full autosomal grid has2897regions.
There are38 TSV exports:12full-grid,24projected,and two RNA-grid identity tables.
These are computational checks, not calibrated genomic intervals or assay resolution.

The unregularized contrast has the positive-count support below. The displayed
RMS range covers all16 technical split deviations from the full, half-count,
median-centred profile, in natural-log units. Positive-count support and finite
display values are distinct; neither certifies biological callability.

| Setting | Query autosomal fragments | Reference autosomal fragments | Positive in both | Split RMS range |
|---|---:|---:|---:|---:|
| trim0_before_MAPQ10 | 1138491 | 2113295 | 2732/2897 | 0.075567–0.087155 |
| trim0_before_MAPQ20 | 1127076 | 2093805 | 2731/2897 | 0.077261–0.093834 |
| trim0_before_MAPQ30 | 334686 | 761028 | 2726/2897 | 0.131316–0.148932 |
| trim0_after_MAPQ10 | 733884 | 1443425 | 2732/2897 | 0.089427–0.105700 |
| trim0_after_MAPQ20 | 725269 | 1429198 | 2731/2897 | 0.092278–0.111030 |
| trim0_after_MAPQ30 | 213204 | 515891 | 2726/2897 | 0.163450–0.184999 |
| trim9_before_MAPQ10 | 2840529 | 4148587 | 2733/2897 | 0.049034–0.057135 |
| trim9_before_MAPQ20 | 2804795 | 4101734 | 2733/2897 | 0.050610–0.058376 |
| trim9_before_MAPQ30 | 1983948 | 3181975 | 2733/2897 | 0.060159–0.070083 |
| trim9_after_MAPQ10 | 1815533 | 2836430 | 2733/2897 | 0.060063–0.072598 |
| trim9_after_MAPQ20 | 1791423 | 2803336 | 2733/2897 | 0.062274–0.078245 |
| trim9_after_MAPQ30 | 1246307 | 2150645 | 2733/2897 | 0.073607–0.093054 |

At the prespecified post-deduplication MAPQ20 contrast, trim0 versus trim9 has
RMS0.103911948, median absolute difference0.047510558 and maximum2.079666062.
Deduplication RMS is0.057050608/0.040057517 for trim0/9. MAPQ20→30 RMS is
0.167266225/0.047176311. These are descriptive processing effects on the same
libraries; all20 fixed comparisons remain in the linked report. No significance,
preferred arm or genotype interpretation is inferred. Half-samples share
libraries and cannot detect a shared amplification/mappability or reference bias.

**Decision: genomic accuracy scoring remains inadmissible.** The frozen gate
requires independently supported reference biology, DNA variation/known-state
information and assay resolution; those measurements are still missing.
No post-outcome RMS cutoff is introduced. NC297 remains an uncertain reference,
and the one exposed patient supplies no independent validation cohort. No RNA
prediction values are parsed and no RNA/DNA correlations, ranks, copy states or
clone labels are computed. Preserve all fixed arms; do not broaden reprocessing
without a specific new discrepancy. More splits or simulations cannot replace
the missing measurements. This chain is now complete, with zero running/pending
jobs and zero failed jobs; the broader scientific goal remains unachieved.

| Stage | Progress | Running / pending / done | Memory and downstream state |
|---|---|---|---|
| Five full RNA/DNA producers | 5/5 ██████████ | 0 / 0 / 5 | Final producer38.9GiB sampled; complete |
| Full native reconstruction | 1/1 ██████████ | 0 / 0 / 1 | Process53068KiB;24 exact comparisons |
| Full DNA sensitivity | 1/1 ██████████ | 0 / 0 / 1 | Process97948KiB;genomic scoring withheld |

[Full result and all20 comparisons](../dev/bian_DNA_profile_sensitivity_20261112/full_23235024/results/REPORT.json),
[all204 checks and every split](../dev/bian_DNA_profile_sensitivity_20261112/full_23235024/results/PROFILE_REPORT.json),
[actual environment](../dev/bian_DNA_profile_sensitivity_20261112/full_23235024/results/ENVIRONMENT.json),
[resources and warnings](../dev/bian_DNA_profile_sensitivity_20261112/full_23235024/resources.txt),
[unaltered warning log](../dev/bian_DNA_profile_sensitivity_20261112/full_23235024/output.txt),
[completion/admission checkpoint](../dev/bian_DNA_profile_sensitivity_20261112/FULL_COMPLETION.json).
Earlier sections retain their submission-time context.

## Why this is the next analysis

The previous turn ruled out the unsupported-reference-RNA correction under its
fixed gate and completed the only second RNA reference's full-library counts.
Its 18 reference/model/baseline outputs are frozen on identical support before
DNA scoring. This is progress, but no new genomic performance improvement.

At the start of this turn, Slurm confirms query trim9 **23233128** complete
0:0 and query trim0 / both reference arms **23233127/23233129/23233130** live.
The independent native reader **23233877** remains pending Dependency. No job
is restarted, moved, cancelled or modified. Source results and historical
decisions remain unchanged. All pending/new HaploDose work uses nslab account
and QOS; unrelated running workloads are untouched.

The earlier [full-measurement contract](BIAN_FULL_MEASUREMENTS_2026-11-12.md)
and [native reconstruction design](BIAN_NATIVE_RECONSTRUCTION_2026-11-12.md)
require an explicit later normalization/support/technical-precision contract
before comparing RNA with DNA. The new analysis implements that requirement
without adding processing arms, a specimen or a biological validation claim.

## Fixed estimand and scientific interpretation

The estimand is query/reference accepted-fragment rate per fixed 1 Mb interval,
under a declared trimming, positional-duplicate and MAPQ setting. Query is
CRC04_LN1_349, reference CRC04_NC_297, distinct cells from one declared patient.
The NC label is an uncertain reference measurement, not a known diploid genotype.
Fragment rates include genome, bisulfite mapping, amplification, GC/mappability
and capture effects. A regional coverage bias or alteration in the reference
can produce a dosage-shaped contrast. Independently corresponding DNA,
replicate libraries and suitable known-state controls would distinguish those
explanations; technical splits of these same reads cannot.

The complete family remains trim0/9 × before/after deduplication × MAPQ10/20/30:
**12 contrasts**, each from two independently reconciled native count tables.
The original post-deduplication MAPQ20 pair is the primary processing contrast;
it is not selected by RNA agreement. All settings, zeros and failures remain.
There is no new RNA-optimized support, segmentation, correction or resolution.

Library totals use **all accepted autosomal fragments** at the same processing
setting, including bins absent from RNA. The unregularized empirical log-rate
ratio exists only where query and reference are positive. Export query-zero,
reference-zero, both-zero and zero-library reasons explicitly. Finite empirical
log ratios do not establish biological callability. Zero libraries are unavailable.

A separate fixed half-count display contrast keeps zero bins finite. Its
full-autosomal-grid median is subtracted for an additional descriptive export.
The half count is not an observed fragment, neutral copy call, justified prior
or calibration assumption. Both zero-count and full-support observations remain
visible. Projection onto either frozen RNA grid preserves the already computed
full-grid totals and centre; no projected support is silently renormalized.

## Exact implementation and falsification

The private [driver](../dev/bian_DNA_profile_sensitivity_20261112/run.py) and
[design](../dev/bian_DNA_profile_sensitivity_20261112/DESIGN.md) freeze **312**
source identities before execution, inheriting both the native-reader and
RNA-lock inventories. No public/default/native package source or test changes.

The full job requires successful native reconstruction and all 24 native/
producer table comparisons, then hashes actual input tables before and after.
It requires the full ordered reference dictionary, exact integer sums,
MAPQ nesting and post-deduplication counts no greater than their predecessors,
including every half-sample. It checks both RNA artifact locks, parses only
their genomic identities and projects without reading RNA prediction values.

All eight deterministic original-qname halves and their complements are retained.
Every bin and library total must conserve half + complement = full exactly.
The 16 split rate profiles use their own library totals. Export every vector,
per-region range, and whole-grid RMS/maximum deviation from the full vector.
These are conditional technical observations, not independent specimens,
biological intervals, uncertainty coverage or a false-split guarantee.

All 20 fixed processing comparisons are retained: six trim pairs, six
deduplication pairs, and eight adjacent-MAPQ pairs. Summaries include finite
support, RMS, median absolute and maximum absolute differences. No stability
cutoff is chosen to declare truth or select a favored processing arm.

Each full and split log contrast is checked against 60-digit Decimal logarithms
of exact integer-rational rates. The arithmetic budget is the previously used
32eps × (1 + maximum effect magnitude), with actual discrepancies reported.
This checks floating-point computation; it is not statistical uncertainty or a
universal libm guarantee. A failed numerical/integrity gate stops the analysis
without relaxing criteria or removing a failed setting.

Even complete success leaves **genomic truth/resolution qualification false**:
an independent biological variation or known-state model is unavailable for this
pair. No RNA/DNA correlation, accuracy rank, reference choice, copy state or clone
label is calculated by this stage. Technical stability cannot repair missing
replication or a shared assay bias. Report those gaps rather than manufacturing
an error threshold from additional compute.

## Executed known-answer qualification

Fixture **23235022 completes 0:0**, nslab/nslab CPU1/2 GB/1 hour,
ne1dc4-011. Slurm elapsed 14 seconds; timed process 10.51 seconds, peak RSS
63,656 KiB. Slurm batch RSS is zero, an unavailable short-job sample, not zero
memory use. All **312** before/after source checks match exactly.

Hand-declared Fraction ratios and independent Decimal logarithms verify unequal
libraries, equal rates, zero counts, zero libraries and median centring.
ID-preserving permutations and projection preserve the original normalization.
All nine invalid cases reject: duplicate grid, reordered grid, missing bin,
broken split conservation, negative count, floating count, nonnested MAPQ,
post-deduplication count inflation and duplicate projection identities.
No biological or statistical calibration follows from these checks.

[Qualification](../dev/bian_DNA_profile_sensitivity_20261112/fixture_23235022/results/QUALIFICATION.json),
[environment](../dev/bian_DNA_profile_sensitivity_20261112/fixture_23235022/results/ENVIRONMENT.json),
[resources](../dev/bian_DNA_profile_sensitivity_20261112/fixture_23235022/resources.txt),
[source inventory](../dev/bian_DNA_profile_sensitivity_20261112/SNAPSHOT.json).

## Submitted full analysis and downstream status

Full job **23235024** requests nslab/nslab CPU1/8 GB/2 hours with
`afterok:23235022:23233877`. It is pending Dependency. Output/result checks are
**not run yet**. It uses the existing isolated read-only Conda/runtime, one
thread and a neutral fresh directory. No installs, arrays, shared writes,
running-job modifications or frozen-output changes.

Exact submissions:

```sh
sbatch --parsable --cpus-per-task=1 --mem=2G --time=01:00:00 --output=/gpfs/commons/home/jameslee/Haplodose/dev/bian_DNA_profile_sensitivity_20261112/slurm-%j.out dev/bian_DNA_profile_sensitivity_20261112/run.sbatch /gpfs/commons/home/jameslee/Haplodose/dev/bian_DNA_profile_sensitivity_20261112 fixture
sbatch --parsable --dependency=afterok:23235022:23233877 --output=/gpfs/commons/home/jameslee/Haplodose/dev/bian_DNA_profile_sensitivity_20261112/slurm-%j.out dev/bian_DNA_profile_sensitivity_20261112/run.sbatch /gpfs/commons/home/jameslee/Haplodose/dev/bian_DNA_profile_sensitivity_20261112 full /gpfs/commons/home/jameslee/Haplodose/dev/bian_DNA_profile_sensitivity_20261112/fixture_23235022
```

Original recovery remains 2/5 in production; the separate 5/5 terminal-pair
research milestone fails the joint error/power requirement and remains
unpromoted. State false splits remain 14/200, simultaneous upper 11.87147%;
matched-generator power 195/200; selected correlated coverage 54/1,000.
These are open scientific failures, not defects repaired by DNA counting.
The fixed competitor family, consumed-data restrictions and all historical
decisions remain. Public pytest/Rust suites are **not rerun** for private
measurement drivers; compute-node fixtures above are executed as reported.

Revisited the [frozen-model failure trace](FROZEN_POPULATION_TRACE_2026-11-11.md)
and [exact-path investigation](EXACT_PATH_MOVES_2026-11-11.md) before considering
another clone intervention. Neither supplies a new better feasible point or
mechanism justifying a change: forcing K2, removing state adjustment and the
tested exact path moves already fail their criteria. No new optimizer sweep,
lowered support or selection-rule rescue is submitted from those old outcomes.

## Final scheduler and preservation checkpoint

Running2, pending2, done3, terminal errors0 in this seven-job chain. At the
last poll the live mappings have elapsed2h04m; neither is near its90h limit.
Memory below is actual sampled current usage for live jobs and timed peak for
the short fixture. Completed producer source checks are257 each; their regional
values remain unopened by this analysis.

| Stage | Progress | Running / pending / done | Memory and downstream state |
|---|---|---|---|
| Four DNA mappings23233127–23233130 | 2/4 ████░░░░ | 2 / 0 / 2 | Current query26.4GiB/reference34.3GiB; native reader waits |
| Native reader23233877 | 0/1 ░░░░░░░░ | 0 / 1 / 0 | Not started; after all four producers |
| Sensitivity fixture23235022 | 1/1 ████████ | 0 / 0 / 1 | 63,656KiB process peak; pass |
| Full sensitivity23235024 | 0/1 ░░░░░░░░ | 0 / 1 / 0 | Not started; after fixture and native reader |

Static AST/shell syntax, ten report links,312 executed before/after source
identities and preservation of all465 preceding ledger values pass. Ledger now
has467entries. These checks preserve computational evidence; no open scientific
failure becomes a pass and no goal-completion claim is made.
