# Full-library measurements toward the genomic profile challenge

Five immutable compute jobs now extend the qualified bounded RNA/DNA processing
to the full prespecified query and reference-candidate libraries. This supplies
new primary measurements for a genomic profile test. No performance result is
available at submission, and the original scientific goal remains unachieved.

## Full chain complete; genomic accuracy scoring withheld

Independent native23233877 now completes0:0,274 source identities and24 exact
count comparisons pass. Sensitivity23235024 completes0:0,312 source identities,
204 independent Decimal checks and all20 processing comparisons pass. The
production/reconstruction/sensitivity chain has no running or pending jobs.
The frozen gate still lacks independently supported reference biology,DNA
variation/known-state information and assay resolution. No RNA predictions
are parsed or RNA/DNA accuracy scores computed by these stages. Both trim arms
and all scientific criteria remain;no performance improvement or biological
validation follows. Earlier checkpoints retain their original historical status.
[Completed measurements, sensitivity and admission decision](BIAN_DNA_PROFILE_SENSITIVITY_2026-11-12.md).

## All five producers complete; independent reconstruction running

All five jobs now complete0:0, with **257 exact before/after source checks per
job**. Final reference trim9 producer23233130 completes in3h31m09s; timed process
elapsed3h31m02s, maximum RSS12808500KiB, versus Slurm batch peak40838552KiB
(38.9GiB). These measure different process scopes and are retained separately.
Its native alignment, quickcheck, positional deduplication, lexical name sorts,
and all six before/after MAPQ10/20/30 tables complete. The five chrM extraction
warnings remain in the unmodified log; no warning is treated as proof of an
autosomal error or silently removed. No running job was restarted or modified.

Independent reader **23233877** starts automatically on nslab/nslab with all274
source-before checks passing. Its full eight role/trim/stage reconstructions
and24 MAPQ comparisons are not complete at this checkpoint. Sensitivity23235024
remains pending Dependency. No regional DNA values are displayed, RNA/DNA score
computed, arm selected or genomic truth admitted. All earlier recovery/error/
power/coverage and competitor criteria remain unchanged. No new batch jobs or
public source changes are needed for this transition.

| Stage | Progress | Running / pending / done | Errors and downstream state |
|---|---|---|---|
| Full RNA and four DNA producers | 5/5 ██████████ | 0 / 0 / 5 | 0; final peak38.9GiB/64GiB |
| Independent native reconstruction23233877 | 0/1 ░░░░░░░░░░ | 1 / 0 / 0 | No reported failure; source-before274pass |
| DNA sensitivity23235024 | 0/1 ░░░░░░░░░░ | 0 / 1 / 0 | Not started; after successful native reconstruction |

[Exact producer/source/resource and scheduler record](../dev/bian_full_measurements_20261112/PRODUCERS_COMPLETE.json).
Earlier checkpoints below retain their original historical status.

## Completed full RNA result and frozen profile continuation

RNA **23233126 completes0:0 in4m32s**. All257 source identities match before
and after. Both complete libraries pass independent preparation, every primary
ID, native all-gene count equality, separately supported candidate-tag counting,
all54,849 annotated genes/zeros, exact annotation and unchanged private index.

| Library | Input pairs | Read eligible | Assigned reads | Assigned genes | Candidate gene/tag pairs |
|---|---:|---:|---:|---:|---:|
|Query CRC04_LN1_349|1,842,166|1,835,846|495,320|8,424|215,028|
|Reference CRC04_NC_426|1,243,474|1,232,480|57,760|4,777|25,874|

Splice-supported assigned reads are157,627 and11,594 respectively. Candidate tags
remain unvalidated molecules. The weak reference, unmapped/ambiguous/nonfeature
categories and all full-run exclusions remain visible; no reference replacement
or claim of raw-count calibration follows. Slurm sampled peak30,743,992KiB and
timed process maximum30,059,236KiB are separate memory summaries.

The unchanged installed RNA estimator subsequently completes as23233604:269
source checks,26 installed-module identities, all-gene Decimal numerical checks
and matrix-axis permutation pass. Six RNA-only predictions/baselines are frozen
at4,636 supported genes/1,630 regions before regional DNA access. All four DNA
jobs remain running; neither DNA support nor genomic accuracy is yet established.

[Full RNA report](../dev/bian_full_measurements_20261112/run_23233126/results/REPORT.json),
[RNA completion/resource record](../dev/bian_full_measurements_20261112/RNA_completion_checkpoint.json),
[frozen prediction protocol and result](BIAN_RNA_PROFILE_LOCK_2026-11-12.md).

Earlier checkpoints below retain their original context, including the point
when all five jobs were running before RNA finished.

## Live execution after the staging gate

Parent23232996 now completes0:0 with all223 before/after source checks, four
archive identities, full paired-QC and native md5sum/gzip checks passing. The
reference candidates supply1,243,474 RNA pairs and13,643,476 DNA pairs. They
remain distinct cells with no certified normality or independent biological n.

All five jobs23233126–23233130 are authoritatively RUNNING. Every257-entry
source-before check passes. All four DNA native known-answer fixtures pass;
RNA's independent read-counter fixture, private-index-before identity, native
tool identities and projected annotation pass. Full count results, source-after
checks and genomic performance remain pending. Initial sampled RNA RSS reaches
30,743,992KiB during native alignment; DNA jobs are still in trimming at that
checkpoint, so their eventual mapping peak is not yet known. All fixed arms
and dependencies remain unchanged.

[Authoritative live checkpoint](../dev/bian_full_measurements_20261112/live_checkpoint.json),
[completed parent](../dev/bian_reference_resume_v2_20261112/completion_checkpoint.json).

Earlier submission-time statements below retain their historical context.

## Scientific question and scope

Can the corresponding query's RNA profile be assessed against independently
constructed DNA evidence at a supported resolution? First establish full gene
read counts and regional DNA fragment counts without selecting a processing arm
using agreement. Then test DNA stability, callability and reference adequacy
before scoring frozen RNA estimators. RNA/DNA observations may disagree because
of processing, amplification, reference composition, transcriptional state or
genomic differences. Additional computation cannot identify the cause by itself.

The exact query remains CRC04_LN1_349, RNA SRR7461687/GSM3241543 and bisulfite
DNA SRR7214590/GSM3153081. The reference candidates remain distinct cells:
RNA CRC04_NC_426/SRR8096032/GSM3443126 and DNA
CRC04_NC_297/SRR7214625/GSM3153116. Their lexical selection is unchanged.
They are not a same-cell reference pair or certified diploid controls.

One patient and three technical cells do not provide independent biological
replication for inference. Bian retains its exposed secondary/descriptive role;
these observations cannot become a new independent method-selection set.
No clone identities, genomic calibration or SOTA claim follows from this stage.

The current public bulk workflow requires at least three normal samples; this
stage supplies one candidate RNA reference. That minimum must not be relaxed,
and repeated reads or cells must not masquerade as independent normal samples.
Any initial one-reference comparison therefore remains a descriptive profile
diagnostic. Validated bulk performance requires a suitable separate panel.

## RNA: full counts with distinct read and tag populations

Process all paired records in original order, requiring exact archive identities
and full-QC pair totals. Keep frozen v1 eligibility for read counts and the
stricter v2 guard for candidate tags. The read branch still includes v1's ACGT
tag screen; it is not a complete census of usable RNA. Candidate gene/tag pairs
remain unvalidated molecules. No molecular calibration is inferred from counts.

The query barcode AACAACCA and reference barcode AGCAGGAA are fixed from the
existing technical-QC reports and matched to the preserved primary RT-primer
table at B14/B27. Verify all eight position marginals in the first 10,000 records.
This supports the candidate extraction layout, not donor identity or diploidy.

Independently reconstruct v1 preparation for every raw record. Reuse the exact
private STAR index, autosomal GTF projection and native binary identities; hash
the full index before and after and verify the projected annotation. Reuse the
same STAR/featureCounts settings without rebuilding the index. Require complete
primary IDs, independent tag counts restricted to v2-eligible records, and exact
native featureCounts equality at all 54,849 genes, retaining zeros. All integer
checks are exact. Emit no normalized profile or DNA score in this measurement job.

## DNA: retain both processing explanations

Each of query/reference runs at fixed front trim 0 and 9. The latter is the
existing adaptation, not established original primer chemistry. Keep cutadapt
4.4, Bismark 0.24.1 nondirectional paired mapping, Bowtie2, positional paired
deduplication and native samtools 1.22.1 lexicographic sorting exactly as in the
qualified pilot. Full cutadapt input totals must match the source QC pair count.

Run the independent known-answer SAM fixture through native conversion/sorting
before real data. Verify all input/reference assets and record native binary
identities/versions and every command. The strict counter retains all rejection
categories and 1-Mb midpoint-fragment bins before/after deduplication at MAPQ
10, 20 and 30; 20 remains the existing primary technical setting. Export the
complete autosomal zero grid and require exact count conservation. Neither
trimming arm can disappear because it agrees poorly with RNA.

These exports are raw fragment measurements. They do not normalize DNA, segment
it, label clones or define qualified DNA copy truth. Independent reconstruction,
resampling and resolution assessment remain prerequisites to interpretation.
When those fail, the corresponding RNA/DNA discrepancy is inconclusive.

## Execution, dependency and current status

All five jobs use account/QOS **nslab**, CPU partition, a 90-hour limit and
`afterok:23232996`. They freeze **257** source/provenance identities. No array,
GPU, interactive QOS, new dependency, shared-installation write or changes to
running jobs. Maximum concurrent request is 40 CPUs and 304 GB, subject to QOS.

| Job | Work | CPUs | RAM | At submission |
|---|---|---:|---:|---|
|23233126|Full RNA, query and reference|8|48 GB|Pending dependency|
|23233127|Full DNA, query trim 0|8|64 GB|Pending dependency|
|23233128|Full DNA, query trim 9|8|64 GB|Pending dependency|
|23233129|Full DNA, reference trim 0|8|64 GB|Pending dependency|
|23233130|Full DNA, reference trim 9|8|64 GB|Pending dependency|

Python static AST and both shell syntax checks pass. Executable full-library
checks are **not run** at submission. Reference continuation 23232996 is live;
its R1 archive identity has passed, while remaining transfer/QC is pending.
No partial file is treated as a complete input. Exact resources, source hashes,
dependency and failures remain in the job records. Earlier failed transfers and
the serialization correction are preserved separately.

## Decisions and unchanged requirements

Wrong identities, fixture failures, native-count disagreement, missing primary
IDs or source/index changes stop dependent work. Empty/failed outputs remain in
the family. Before any RNA/DNA scoring, lock estimator/baseline family, alignment,
normalization, DNA callability and resolution checks. Resampling cannot validate
unmodelled amplification bias or create missing independent normal specimens.

The frozen wheel, public defaults, historical outputs and all prior claim values
remain unchanged. Recovery remains 2/5, continuous-state false splits 14/200,
narrow-generator power 195/200 and selected correlated coverage 54/1,000. The
full fixed comparator and biological-validation requirements remain. Neither a
single successful profile nor additional compute establishes superiority.

[Locked measurement contract](../dev/bian_full_measurements_20261112/DESIGN.md),
[RNA implementation](../dev/bian_full_measurements_20261112/rna.py),
[DNA implementation](../dev/bian_full_measurements_20261112/dna.py),
[exact DNA changes from qualified pilot](../dev/bian_full_measurements_20261112/DNA_DELTA_RECORD.json),
[RNA resource header](../dev/bian_full_measurements_20261112/rna.sbatch),
[DNA resource header](../dev/bian_full_measurements_20261112/dna.sbatch),
[actual submissions](../dev/bian_full_measurements_20261112/submissions.json),
[source snapshot](../dev/bian_full_measurements_20261112/SNAPSHOT.json),
[live reference chain and preserved failures](BIAN_REFERENCE_CANDIDATES_2026-11-12.md),
[qualified read/tag correction](BIAN_READ_TAG_EXPORT_2026-11-12.md),
[qualified DNA pilot](BIAN_DNA_READ_QC_2026-11-12.md).
