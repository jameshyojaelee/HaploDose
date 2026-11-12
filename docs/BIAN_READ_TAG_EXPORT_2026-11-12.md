# Preserve read evidence separately from candidate molecular tags

The new private schema-3 export restores **46 gene-assigned reads**, increasing
the first-block total from **26,780 to 26,826**. Every gene's candidate-tag count
is unchanged, totalling **22,364**. This is a demonstrated measurement-accounting
correction; it establishes no CNA accuracy, clone recovery or calibrated
uncertainty benefit. Frozen exports and public package source are unchanged.

## Cause and smallest justified change

The v2 primer-overlap guard appropriately excludes questionable molecular tags.
Using that same eligibility for read counts also discarded aligned transcript
evidence. The prespecified matched ablation found 46 assigned reads, including
28 splice-supported reads, in exactly 9,000 rejected-tag records. Most of that
stratum remains unmapped; no important biological effect size is inferred.

`measurements.prepare` now returns two explicit evidence populations: frozen v1
eligibility for reads and frozen v2 eligibility for candidate tags. All other v1
filters, including its ACGT-tag screen, remain. This is not comprehensive
recovery of all potentially useful RNA reads. Tag text in read names provides
opaque read identity only. Exact candidate gene/tag pairs remain distinct from
validated molecules and cannot establish an independent-count sampling law.

## Independent falsification and result

Job **23232893 completed 0:0**, account/QOS nslab, CPU partition, one CPU,
2 GB requested and one-hour limit. It executed for 18 seconds; timed peak RSS
is 252,708 KiB. Assigned reads now span 4,292 genes. Before/after
checks agree for all **217 frozen identities**, including transitive source,
producer reports and the preserved failed launch. The installed isolated
environment and captured runtime loader are used with `python -I`.

All unchanged preprocessing fixtures pass. For every one of the first 100,000
raw paired records, independently written v1 preparation equals the read arm,
and frozen v2 preparation equals the tag arm. Parent and rejected-stratum FASTQs
are disjoint and reconstruct the new **99,665-read** eligibility set exactly;
the parent alone exactly reconstructs the unchanged **90,665-tag** set.

An independent parser of both native assigned SAMs verifies every primary ID
exactly once. Its counts equal the sum of the two featureCounts tables at every
one of **54,849** genes, retaining annotated zeros. Every candidate-tag gene
count is identical to the parent. Integer equalities require exact agreement;
no tolerance or scientific threshold changes. Original raw inputs are rehashed
on the compute node, and transitive FASTQ/SAM/table identities are recorded.
No remapping, index rebuild, feature selection or normalization is performed.

## Preserved failed attempt

Job **23232887 failed 1:0 before analysis**, in zero recorded seconds. The local
preparation command used unavailable `python`, and the submission shell did not
stop before submitting without `source.sha256`. Its launcher correctly refused
to proceed. The correction uses installed `python3` and `set -euo pipefail`
before snapshot creation and submission. The same analytical code, eligibility,
resources and criteria were then submitted as 23232893. Both logs and the
failure record remain; no failed artifact was overwritten or counted as a pass.

## Interpretation and remaining evidence

One technical cell supplies no independent biological replication (biological
inference n=0). RNA-origin alternatives, tag chemistry, normal-reference
adequacy, corresponding DNA callability and genomic membership remain unresolved.
The original recovery 2/5, state false splits 14/200, selected correlated coverage
54/1,000 and full comparator criteria remain unchanged. Bian retains its exposed,
secondary/descriptive role. This export cannot promote a genomic-validation
claim or turn integer counts into calibrated uncertainty.

[Locked design](../dev/bian_read_tag_export_20261112/DESIGN.md),
[eligibility implementation](../dev/bian_read_tag_export_20261112/measurements.py),
[independent reconstruction](../dev/bian_read_tag_export_20261112/run.py),
[actual resource header](../dev/bian_read_tag_export_20261112/run.sbatch),
[completed report](../dev/bian_read_tag_export_20261112/run_23232893/results/REPORT.json),
[schema and support declarations](../dev/bian_read_tag_export_20261112/run_23232893/results/SCHEMA.json),
[transitive artifact hashes](../dev/bian_read_tag_export_20261112/run_23232893/results/ARTIFACT_IDENTITIES.json),
[failed-launch record](../dev/bian_read_tag_export_20261112/failed_submission_23232887.json),
[preceding native ablation](BIAN_TAG_READ_ABLATION_2026-11-12.md).
