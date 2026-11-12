# Native DNA reconstruction and the November10 review

The review's concrete source corrections are already implemented and remain
consistent with the protected source snapshot. Its scientific objections remain
open under their original criteria. This follow-up adds an independently written
native fragment check before interpreting the newly mapped DNA; it does not
promote a public method or claim an accuracy improvement.

## Completed full reconstruction

Full job **23233877 completes0:0**, nslab/nslab CPU2/8GiB/90h,
ne1dc4-023. Slurm elapsed19m35s; timed process19m28.83s, peak53068KiB;
Slurm sampled peak35864KiB. All274 source identities match before/after.
All eight query/reference × trim0/9 × before/after reconstructions and their
24 MAPQ10/20/30 comparisons agree exactly with the producer, including
complete zero grids, nonautosomes, coordinate boundaries and total support.
All eight original-qname hash halves and complements conserve counts exactly.

Dependent sensitivity23235024 subsequently completes0:0 with312 exact source
identities and all204 independent Decimal checks. No mismatch requires a
source correction or broader reprocessing. These results establish the declared
computational measurement contract. They do not establish biological independence
of fragments, an independent DNA assay, copy-state truth or assay resolution.
RNA predictions remain unopened by these stages; genomic accuracy admission
remains false under the previously frozen requirements.

[Full native report](../dev/bian_native_reconstruction_20261112/run_23233877/results/REPORT.json),
[native resources](../dev/bian_native_reconstruction_20261112/run_23233877/resources.txt),
[completed sensitivity and scientific decision](BIAN_DNA_PROFILE_SENSITIVITY_2026-11-12.md),
[complete source/numerical/scheduler checkpoint](../dev/bian_DNA_profile_sensitivity_20261112/FULL_COMPLETION.json).
Earlier checkpoints below retain their historical status.

The [17-finding response](SCIENTIFIC_REVIEW_RESPONSE_2026-11-10.md) separates
implemented corrections from unestablished claims. Production recovery remains
2/5, continuous-state false splits14/200 with simultaneous upper11.87147%, and
selected correlated coverage54/1000. Technical agreement or completed jobs
cannot turn those failures into passes.

## Implemented route and exact contract

New private [reader](../dev/bian_native_reconstruction_20261112/run.py),
[frozen design](../dev/bian_native_reconstruction_20261112/DESIGN.md), and
[launcher](../dev/bian_native_reconstruction_20261112/run.sbatch).
It imports none of the producer's CIGAR/classification functions. A separate
paired-record scanner selects intact proper pairs with known MAPQ and reciprocal
mate metadata. Samtools1.22.1 encodes the selected SAM; BEDTools2.31.0 supplies
native reference-consuming mate coordinates. Exact integer midpoint counting
must reconstruct every producer bin atMAPQ10/20/30, before/after deduplication,
query/reference and both trim0/9 arms. Zero bins and nonautosomes remain in the
comparison. A mismatch or unavailable arm stops the reader rather than dropping
an output from the fixed family.

The full reader also writes eight deterministic original-qname half-samples
and complements. Every bin must conserve counts exactly. These assess only
conditional technical variability. They cannot measure amplification bias,
independent-library variation or biological error. Original qnames share the
same split across processing arms. Reads are not biological replicates.

This route shares the alignment, reference and acceptance definition with the
producer. It independently tests coordinate arithmetic and implements support
eligibility separately; it is not an independent DNA assay or universal SAM
validator. Implementations agreeing could still share a scientifically wrong
measurement definition.

Primary documentation confirms BEDPE's paired-name grouping, coordinate
construction and minimum-mate MAPQ behavior:
[BEDTools2.31](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html),
[samtools](https://www.htslib.org/doc/samtools.html).
Actual compute-node executable versions and SHA256 are recorded below.
No production dependency is added or installed.

## Executed qualification and pending full run

Qualification **23233873 completes0:0**, CPU node
`ne1dc6-008.nygenome.org`, nslab account/QOS,1CPU2GB1h requested. Slurm elapsed
6s; timed scientific process1.83s, maximum RSS53,104KiB. Slurm sampled batch RSS
is0, an unavailable short-job sample, not evidence of zero memory consumption.
All **274** source identities match before/after.

All15 original hand-counted template fixtures remain, including unknown MAPQ,
broken mate correspondence, genomic splicing, duplicates, excess primary records
and an auxiliary alignment. Native coordinates of all seven eligibleMAPQ10
templates agree exactly with their hand-declared intervals. AtMAPQ20/30 there
are four autosomal fragments (threebin0,onebin1) and onechrX fragment. AtMAPQ10
there are six autosomal fragments (fivebin0,onebin1) and onechrX fragment.
All three thresholds and explicit empty autosomal bins pass. No numerical
tolerance is used: the expected objects are exact integer counts.

Evidence:
[qualification](../dev/bian_native_reconstruction_20261112/run_23233873/results/QUALIFICATION.json),
[native tools](../dev/bian_native_reconstruction_20261112/run_23233873/results/TOOLS.json),
[resource record](../dev/bian_native_reconstruction_20261112/run_23233873/resources.txt),
[immutable source inventory](../dev/bian_native_reconstruction_20261112/SNAPSHOT.json).

Full reader **23233877** is submitted with
`afterok:23233873:23233127:23233128:23233129:23233130`.
It requests nslab/nslab CPU2CPU8GB90h, no array. At this checkpoint it is pending
for those dependencies, while all four producers remain running. There is no
full reconstruction result yet, and no RNA/DNA correlation has been computed.
Producer source, native settings and original running allocations are unchanged.

Exact submission commands:

```sh
sbatch --parsable --cpus-per-task=1 --mem=2G --time=01:00:00 --output=/gpfs/commons/home/jameslee/Haplodose/dev/bian_native_reconstruction_20261112/slurm-%j.out dev/bian_native_reconstruction_20261112/run.sbatch /gpfs/commons/home/jameslee/Haplodose/dev/bian_native_reconstruction_20261112 qualification
sbatch --parsable --dependency=afterok:23233873:23233127:23233128:23233129:23233130 --output=/gpfs/commons/home/jameslee/Haplodose/dev/bian_native_reconstruction_20261112/slurm-%j.out dev/bian_native_reconstruction_20261112/run.sbatch /gpfs/commons/home/jameslee/Haplodose/dev/bian_native_reconstruction_20261112 full /gpfs/commons/home/jameslee/Haplodose/dev/bian_native_reconstruction_20261112/run_23233873
```

## Scientific interpretation and next decision

Estimand here: accepted DNA fragments per declared genomic bin under specified
processing, not absolute copy states. The count-generating mechanism includes
bisulfite mapping, library amplification, positional duplicate handling, GC and
mappability. The strongest alternative to dosage is region-dependent coverage
bias combined with an uncertain reference cell. An NC title does not establish
diploidy or single-CNA identity. Technical splitting cannot distinguish these
mechanisms; suitable independent known-state or replicated DNA is still needed.

Raw and normalized regional coverage must remain distinct from validated truth.
The design fixes any later descriptive contrast's library totals to all accepted
autosomal fragments at the same setting, but deliberately admits no genomic
scoring at this stage. A later scorer requires a separately frozen exact
normalization, support and uncertainty contract. No post-DNA RNA estimator,
reference, processing-arm or resolution selection is permitted. Grid1Mb is not
demonstrated assay resolution. Both processing arms remain visible, including
failures. Agreement on this exposed patient is regression/descriptive evidence.

No new public package changes are justified merely by the review's unresolved
biological limitations. All historical decisions and frozen outputs remain.
Public pytest/Rust suites: **not rerun** for this private reader; compute-node
native known-answer qualification executed as above. Full-reader checks remain
pending. No clone-validation, uncertainty, performance or completion claim.
