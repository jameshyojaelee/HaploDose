# Frozen RNA predictions before the corresponding DNA assessment

The installed HaploDose estimator and matched naive baselines now produce locked
RNA profiles from the new complete read counts. Job **23233604 completes0:0**;
all269 before/after source identities, all26 installed module hashes, exact
count/coordinate checks, independent60-digit numerical contrasts and reversed
matrix axes pass. No regional DNA measurement was opened or scored. This is
execution toward the genomic profile challenge, not a demonstrated accuracy gain.

## Scientific contract

The estimand is median-centred natural-log relative RNA abundance in supported
genomic regions against the specified reference. Observations are assigned read
counts from the qualified complete libraries, not independently validated UMIs.
The working smoother uses conditional masked-gene loss. RNA state, reference
composition, RNA quality and genomic dosage remain competing explanations. A
smaller fitted RNA loss does not independently validate the selected estimator
or identify a genomic alteration.

Query CRC04_LN1_349 supplies495,320 assigned reads; reference CRC04_NC_426 supplies
57,760. The reference remains one uncertain candidate cell, with materially less
gene support. It is not replaced based on yield. Candidate tag counts never
substitute for read counts or acquire molecule independence. The public bulk
minimum of three normal samples remains unchanged; this is a descriptive
single-query profile, not a calibrated bulk analysis or a clone experiment.

Biological inference n=0, one exposed patient. The DNA reference candidate is a
different cell from the RNA reference. Bian retains its secondary/regression-only
role. Original recovery/error/power/coverage, comparator and prospective-validation
requirements remain. No scientific or public/default criterion changes.

## Fixed features and fair baselines

Gene records come from the exact Ensembl GRCh37r75 GTF. Convert intervals once
to0-based half-open coordinates and verify hg19 chromosome bounds. The existing
UCSC hg19 cytoband archive matches its source-manifest MD5; each centromeric gap
is the union of two adjacent acen bands. Retain positive-reference autosomal
genes with explicit symbols, excluding the unchanged public ribosomal category
and genes/grid bins overlapping the centromere. This support rule matches the
usual1CPM minimum for this57,760-read reference and uses no query response or DNA.

**4,636 genes in1,630 occupied1Mb regions** remain. Exclusions:50,072 genes
without reference reads,107 ribosomal genes and34 genes whose gene/grid support
overlaps the centromere. Selected totals are352,213 query reads and51,085 reference
reads. Missing reference support remains missing; no region is filled to improve
coverage. The1Mb grid is not validated biological resolution. Blocks split at
chromosome/arm boundaries and every missing physical bin; terminal bins are
clipped to the exact chromosome length.

All three methods use identical selected genes, reference normalization and
public half-count convention: naive within-bin gene mean, naive within-bin
gene median and the unchanged partition-average RNA stack. Export each both
as returned and under the identical width5 within-block median operator.
These six locked predictions distinguish adaptive-estimator behaviour from a
common smoothing effect. The one query cannot supply a target-excluded
pseudobulk population; these are single-cell ratio baselines, not a replacement
of the original pseudobulk benchmark.

The existing noneditable0.2.0 wheel runs
`stacked_count_profile(...,counts_kind='raw',gene_segments=True,average_segments=True,solver='primal')`.
Every installed package Python module matches the earlier software freeze.
Its partition-average expert is selected by its existing RNA criterion; that
selection is conditional fitting, not a passed genomic-performance criterion.
No model, prior, variance, native dependency or public source is changed.

## Independent computational checks

Every feature identity and integer subtotal is conserved. All gene contrasts
are independently reconstructed with Decimal60-digit logs of exact rational
pseudocount ratios, followed by the same median centre. Maximum absolute error
is **8.881784197001252e-16**, below the input-dependent operation budget
**6.92689243935283e-14**. The32eps magnitude expression accounts for the rate,
offset, ratio, log and centring arithmetic; it is a numerical check under ordinary
floating-point assumptions, not a universal libm theorem or statistical interval.

Reversing both matrix axes, preserving IDs and the fixed feature order, gives
bitwise-identical gene effects and profiles. The common median operator is
independently implemented by sorted bounded windows. Prediction, feature and
mask hashes are locked together with original count/GTF/cytoband identities and
resolved method diagnostics. No future DNA-file path is read by this script.

Job23233604 uses nslab account/QOS, CPU1CPU8GB1h and a neutral cwd in the existing
read-only isolated environment. Slurm elapsed1m11s; timed process1m07.91s.
Timed peak RSS199,008KiB, Slurm sampled155,580KiB. Two SciPy OptimizeWarnings
about the `threads` option passed through to HiGHS are retained. They are not
silenced or treated as proof of solver accuracy. All269 source identities pass
before and after; the complete public pytest/Rust suites are not rerun.

## Next decision and falsification

Four full-library DNA jobs remain live, retaining query/reference × trim0/9,
MAPQ10/20/30 and before/after-deduplication measurements. Before any RNA/DNA
score, independently reconstruct DNA fragments, assess processing stability,
resampling/callability and the supported resolution. A stable technical coverage
profile alone cannot exclude amplification bias or certify normality.

If suitable DNA support fails, report the genomic question inconclusive. If a
profile is constant or missing, preserve its unavailable status under the
unchanged common-support rule. A difference confined to unmatched smoothing
cannot establish the proposed mechanism. Any one-patient agreement remains
descriptive; no confidence coverage, clone membership or SOTA claim follows.
The full goal remains active and unachieved.

[Locked design](../dev/bian_rna_profile_lock_20261112/DESIGN.md),
[executed RNA-only implementation](../dev/bian_rna_profile_lock_20261112/run.py),
[native resource header](../dev/bian_rna_profile_lock_20261112/run.sbatch),
[actual job record](../dev/bian_rna_profile_lock_20261112/submission.json),
[frozen predictions and diagnostics](../dev/bian_rna_profile_lock_20261112/run_23233604/results/RNA_LOCKED.json),
[prediction table](../dev/bian_rna_profile_lock_20261112/run_23233604/results/PREDICTIONS.tsv),
[all numerical checks](../dev/bian_rna_profile_lock_20261112/run_23233604/results/NUMERICAL_CHECKS.json),
[coverage/mask](../dev/bian_rna_profile_lock_20261112/run_23233604/results/MASK.json),
[completed source/resource checkpoint](../dev/bian_rna_profile_lock_20261112/completion_checkpoint.json),
[full source measurements and DNA jobs](BIAN_FULL_MEASUREMENTS_2026-11-12.md).
