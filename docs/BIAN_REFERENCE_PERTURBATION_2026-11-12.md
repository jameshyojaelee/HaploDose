# Bian: preserve the weak reference and test the only second RNA reference

The second CRC04 NC-labelled RNA library now has independently reconciled
full-library read counts. This extends a fixed RNA reference-sensitivity
comparison; it supplies no new patient, diploid standard, genomic clone truth
or validated uncertainty. The original reference, frozen outputs and scientific
criteria remain intact.

## Why this addition, and exactly which material

The historical RNA reference CRC04_NC_426 has 57,760 assigned reads and 4,777
genes. The frozen primary ENA metadata contains exactly two NC-labelled RNA
cells for CRC04: 426 and 429. The new library is the only second cell, rather
than a reference chosen by RNA fit or DNA agreement. No third cell or patient
is added to satisfy the bulk minimum of three normals.

The selected row is `scTrioSeq2Rna_CRC04_NC_429 [single cell RNA-seq]`,
GSM3443127 / SRR8096033, RNA-Seq / TRANSCRIPTOMIC, paired HiSeq4000 under the
declared scTrioSeq2 protocol. There is no corresponding NC same-cell DNA pair
in the previous 609-pair crosswalk. Metadata correspondence does not establish
donor concordance, diploidy or lineage matching. Both references remain distinct
from DNA reference CRC04_NC_297.

[Frozen design](../dev/bian_reference_perturbation_20261112/DESIGN.md),
[archive identities](../dev/bian_reference_perturbation_20261112/INPUTS.json),
[executed code](../dev/bian_reference_perturbation_20261112/run.py).

## Attempt 1: archive transfer and complete paired QC

Job **23234650 completes 0:0**, nslab/nslab I/O1 CPU/2 GB/2 hours,
ne1di6-002, Slurm elapsed 1m01s; timed 1m00.35s, process peak 41,788 KiB,
scheduler sampled peak 27,548 KiB. The I/O
partition is used for transfer and hashing. Both provider bytes and MD5, native
`md5sum`, `gzip -t`, full paired-record QC and malformed FASTQ fixtures pass.
All **277** source identities match before and after.

| File | Provider bytes | Provider MD5 |
|---|---:|---|
| SRR8096033_1.fastq.gz | 69,784,474 | ff7cf384dd6e60ba13acf5d35ecc519e |
| SRR8096033_2.fastq.gz | 74,950,743 | 2772b147731f423c208e53b9db9fb188 |

Total archive bytes are 144,735,217; complete paired records are 976,074.
The exact HTTPS archive paths and SHA256 values are saved with the technical
input checks. This is measurement integrity, not a molecular or genomic claim.

[Stage result](../dev/bian_reference_perturbation_20261112/stage_23234650/results/REPORT.json),
[full read QC](../dev/bian_reference_perturbation_20261112/stage_23234650/results/RNA_QC.json),
[native checks](../dev/bian_reference_perturbation_20261112/stage_23234650/results/NATIVE_CHECKS.json).

## Attempt 2: unchanged full-library native and independent counting

Job **23234654 completes 0:0**, afterok 23234650, nslab/nslab CPU8/48 GB/
90 hours, ne1dc5-011. Slurm elapsed 4m47s; timed process 4m42.43s, peak
29,975,000 KiB; scheduler sampled peak 30,289,676 KiB. These are different
process/scheduler measurements. All **277** before/after identities match.

The code reuses the qualified full-library counter unchanged: private hg19
STAR 2.7.11b index, Ensembl75 autosomal exon projection, featureCounts 2.1.1,
samtools 1.22.1 and the existing read/tag preparation. Index identities match
before and after; no index rebuild, mapping tuning, installation or shared write.

The first eight R2 positions are pure across the first 10,000 technical records
and match exactly one primary Barcode-RT entry: **ATTGAGGA, B30**. Pure-unique
fixtures pass; mixed, unknown-base, absent and duplicate matches are rejected.
That supports this technical layout, not specimen or molecule identity.

| Full-library quantity | Count |
|---|---:|
| Input pairs | 976,074 |
| Eligible read observations | 970,173 |
| Candidate-tag-eligible observations | 915,138 |
| Assigned reads | 94,056 |
| Genes with assigned reads | 4,044 |
| Candidate gene/tag pairs | 30,384 |
| Assigned reads with splice CIGAR | 32,459 |
| Annotated genes retained, including zeros | 54,849 |

Independent all-record preparation, exact primary IDs, every native gene-read
count and candidate-tag count agree. The read-counter rejection fixtures pass.
All status categories are retained: 94,056 assigned, 770,877 unmapped, 36,120
nonunique, 66,487 secondary/supplementary, 57,880 no feature and 11,240 ambiguous.
Candidate tags remain unvalidated molecules; integral read counts do not justify
raw-UMI calibration or independent Poisson/NB sampling assumptions.

[Complete count result and executed native commands](../dev/bian_reference_perturbation_20261112/counts_23234654/results/REPORT.json),
[barcode](../dev/bian_reference_perturbation_20261112/counts_23234654/results/BARCODE.json),
[all gene counts](../dev/bian_reference_perturbation_20261112/counts_23234654/results/reference429/GENE_COUNTS.tsv),
[input archive hashes](../dev/bian_reference_perturbation_20261112/counts_23234654/results/reference429/INPUT_CHECKS.json),
[resources](../dev/bian_reference_perturbation_20261112/counts_23234654/resources.txt).

## Attempt 3: fixed reference comparison, submitted before genomic scores

The comparison was specified in the perturbation design before new counts.
After counts qualify, its exact executable/source inventory freezes **298**
identities. Job **23234786 completes 0:0** on nslab/nslab CPU1/8 GB/2 hours,
original profile host ne1dc4-039. Slurm elapsed 1m46s; timed process 1m42.99s,
peak RSS 155,160 KiB; scheduler sampled peak 122,564 KiB. All **298** before/
after identities match. Seven existing SciPy/HiGHS option-forwarding warnings
are retained; no options are changed or warnings suppressed to obtain a pass.

Query stays CRC04_LN1_349. Use the original 4,636-feature mask intersected with
positive 429 reads, without selecting a new coverage cutoff. All three reference
constructions use that same ordered support, hg19 coordinates, physical 1 Mb
bins and block splits at physical/arm/chromosome gaps. Drop empty bins explicitly.
The historical 426 original-support output is preserved and replayed exactly;
the reprocessed 426 arm on common support is a separate matched counterfactual.

Fixed references: 426 alone, 429 alone and their raw-read-weighted pool. For
each, the frozen installed estimator, naive gene mean and naive median are
exported both natively and under the same within-block median5: **18 outputs**.
The single query cannot supply a target-excluded pseudobulk population. No
baseline is renamed to suggest otherwise; no best reference is selected.

Required computational gates include all 26 installed module hashes, bitwise
historical replay, exact ID-preserving axis permutations, and independent
60-digit Decimal gene effects against the original input-magnitude operation
budget. All profiles lock before regional DNA scoring. Failed integrity or
numerical checks reject the outputs; they are not repaired by relaxing tolerances.

All required gates pass. Historical 426 original-support predictions and gene
effects reproduce bitwise. Every reference is invariant to ID-preserving matrix
axis reversal. The common mask contains **2,398 genes / 1,166 regions**, excluding
2,238 original genes without positive 429 support. This support loss remains
explicit and is not evidence that excluded regions are genomically neutral.

| Fixed reference | Query selected reads | Reference selected reads | Maximum gene-effect discrepancy | Arithmetic budget |
|---|---:|---:|---:|---:|
| 426 | 312,789 | 37,615 | 1.78e-15 | 7.06e-14 |
| 429 | 312,789 | 60,997 | 8.88e-16 | 6.53e-14 |
| Read-weighted 426+429 | 312,789 | 98,612 | 8.88e-16 | 6.78e-14 |

The budgets describe numerical arithmetic conditional on the fixed convention,
not inferential coverage. All **18 RNA-only profiles**, features, decisions and
input/output SHA256 values lock before DNA scores. No reference is selected from
these outputs. Neither improved biological accuracy nor reduced uncertainty has
been tested by this successful execution.

[Exact computational checks](../dev/bian_reference_profiles_20261112/run_23234786/results/CHECKS.json),
[locked inputs, outputs and model decisions](../dev/bian_reference_profiles_20261112/run_23234786/results/RNA_LOCKED.json),
[all 18 predictions](../dev/bian_reference_profiles_20261112/run_23234786/results/PREDICTIONS.tsv),
[common support](../dev/bian_reference_profiles_20261112/run_23234786/results/COMMON_FEATURES.tsv),
[environment](../dev/bian_reference_profiles_20261112/run_23234786/results/ENVIRONMENT.json),
[resources](../dev/bian_reference_profiles_20261112/run_23234786/resources.txt).

The estimand is median-centred relative RNA on common support against a declared
reference. Stable effects may reflect dosage or stable transcription; reference
sensitivity may reflect noise or composition. Independent corresponding DNA,
replicated RNA and reference biology are needed to distinguish them. No model
selection, calibrated interval, DNA state or clone claim follows from this stage.

[Exact comparison design](../dev/bian_reference_profiles_20261112/DESIGN.md),
[implementation](../dev/bian_reference_profiles_20261112/run.py),
[submission header](../dev/bian_reference_profiles_20261112/run.sbatch).

## Unchanged scientific decisions

Bian remains exposed secondary/descriptive evidence. Biological specimens added:
zero. Bulk minimum normals: three. The four existing DNA mappings and their
dependent native reader remain untouched; no DNA region values or RNA/DNA
scores are opened here. Resolution-matched, resampling-stable DNA still precedes
genomic interpretation, and stability alone cannot certify biological truth.

Original recovery 2/5, state-null false splits 14/200, selected correlated
coverage 54/1,000, frozen comparator families and all original criteria remain.
No package/default/native/test change or frozen-artifact overwrite is made.

Final operational checkpoint: query trim9 **23233128 completes 0:0** in
1h47m59s, scheduler sampled peak 29,369,616 KiB, all257 before/after source
identities exact. Jobs23233127/23233129/23233130 remain running;23233877 is
pending Dependency. No regional DNA values are read to make this status report.

Static/preservation checks: four driver ASTs, five shell launchers, 32 local
report links and1,096 completed compute source checks pass. All456 preceding
ledger values remain exact;465 entries now record both positive and negative
results. [Checkpoint](../dev/bian_reference_profiles_20261112/preservation_checkpoint.json).
