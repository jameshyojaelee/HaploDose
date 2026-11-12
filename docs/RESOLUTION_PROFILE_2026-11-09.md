# Fixed physical resolution does not beat the retained bulk predictor

Both implemented10Mb-grid RNA candidates fail the bulk performance comparison.
Predictive smoothing on this grid reduces mean correlation from .680245 to
.663768 and loss recall from .311493 to .256273. The unsmoothed grid is worse
for correlation and precision. Reject both without retuning the grid. No default,
frozen result, comparator, scientific criterion or production file changes.

## Intervention and scientific contract

The question was whether fixed gene-count bins, with unequal physical spans,
explain the gap to strong approximately10Mb DNA-profile comparators. The estimand
is descriptive regional relative RNA abundance, using the existing natural-log
gene contrasts, references and RNA-selected gene universe. This is conditional
median prediction, not a calibrated NB count model, DNA state or clone estimate.
Regional transcriptional programmes and reference composition remain alternatives
to genomic dosage. Corresponding independent DNA at the declared resolution and
new biological replication would establish an external performance benefit.

`resolution_profile_panel.py` assigns each existing RNA bin by its midpoint to
fixed10Mb hg38 intervals within the same arm/gap blocks. Genes remain in their
original bins at boundaries. Missing physical intervals and declared gaps stop
smoothing. The rule creates284 coarse intervals from1131 fine intervals.
No DNA values, DNA segments,
purity, labels or accuracy-selected resolution enter grid construction.

Two candidates use either the coarse median alone or the unchanged conditional
masked-gene median selector on those intervals. The RNA prediction NPZ, full
selection records and hashes are saved before DNA opens. Independent hand-built
boundary/gap and noiseless fixtures pass; every installed predictive RNA profile
reproduces exactly for all68 samples before scoring. Gene-anchored DNA mapping
still uses hg19 while RNA coordinates are hg38. Nominal scale matching does not
establish exact physical correspondence, and uncertain published truth remains
uncertain.

## Complete fixed comparison

The executed archive records1131 fine RNA bins and284 coarse intervals. All
twelve original comparators plus fine median1, median5 and installed predictive
RNA are retained:17 total methods including both candidates. Existing common
support,24/24 tumour split, thresholds and missing-patient rules remain unchanged.
There are48 existing tumour patients and20 shared normal reference donors;
the68 RNA samples are not68 independent tumour patients. All data were exposed
previously and remain regression-only.

| Method | Mean correlation (n=23) | Loss recall (n=13) | Loss precision (defined n) | Loss FPR (n=23) |
|---|---:|---:|---:|---:|
| Retained gene predictive |.680245|.311493|.284772 (16)|.042193|
|10Mb median |.597304|.338736|.151989 (21)|.049060|
|10Mb predictive |.663768|.256273|.339259 (13)|.042233|
| Strong naive own partition |.730261|.399679|.358571 (14)|.041244|
| Existing precision-weighted HD partition |.733069|.577106|.458667 (15)|.046896|

Against retained predictive RNA,10Mb median correlation changes -.082942
(5wins/18losses) and recall +.027244. Its worse precision/FPR prevents interpreting
the recall rise as a shared benefit.10Mb predictive changes correlation -.016477
(10wins/13losses), recall -.055220 (0wins/4losses/9ties), and FPR +.0000394.
Its apparent precision rise has fewer defined patients; on the13 paired defined
patients the change is +.027232. This does not compensate for lost recovery.
All no-call/undefined outcomes and every comparator-specific denominator remain
in `patient_scores.tsv`, `summary.tsv` and `report.json`.

The new models also lose against strong naive segmentation, by -.132958 and
-.066493 mean correlation, respectively. These outcomes reject physical
coarsening as the proposed sufficient performance correction, without proving
that arbitrary segmentation or resolution is biologically identified. No
uncertainty-calibration or independent-superiority claim follows.

## Execution and preservation

Job23144495 completes0:0 in42s on nslab CPU, requesting4CPUs/16GB/2h. The task
retains dirty baseline, immutable source_v1, exact submission, installed wheel,
neutral-directory imports, process/scheduler resource records and full output
hashes. Sampling reports102524KiB batch RSS; process peak is separately retained
and can exceed the scheduler sample. The current noneditable artifact's prior
412 Python/4 Rust passes and1 GPU-only skip are separate evidence; full production
and GPU suites were not rerun for this research-only file.

[Frozen design](../dev/resolution_profile_20261109/DESIGN.md),
[complete comparisons](../dev/resolution_profile_20261109/bulk_v1/report.json),
[sample-level scores](../dev/resolution_profile_20261109/bulk_v1/patient_scores.tsv).
