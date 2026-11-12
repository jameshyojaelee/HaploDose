# Andor DNA fetchers misplaced reads at exact bin boundaries

Both Andor DNA fetchers now convert one-based SAM alignment starts to zero-based
500 kb bins with `floor((POS-1)/500000)`. Previously they used `floor(POS/500000)`,
placing a start exactly divisible by 500000 in the following bin. The independent
audit confirms ten such filtered reads in prespecified bounded archive queries.
Historical counts, DNA labels and RNA scores remain unchanged. The effect on
genome-wide DNA profiles or clone reconstruction is not established here.

## Scientific contract and measurement chain

The corrected quantity is a filtered fragment-start count per read group and
zero-based half-open genomic bin. It is a DNA read-depth measurement, not an RNA
model or a direct copy-number observation. The original read filter retains
proper-pair read1, removes flags 0xF04 and requires MAPQ at least30. Both fetchers
retain that filter and their read-group handling. Counts are conserved by the
coordinate mapping; exact-boundary starts move one bin to the left.

The locked development lines remain SNU-601, MKN-45 and NUGC-4. The historical
source trace independently reconstructed16 inferred DNA groups in SNU-601, with
median236973 callable fragments, versus one inferred group at median3220 and5665
in MKN-45 and NUGC-4. Those are inferred measurements: the latter two are not
certified single-CNA populations. Separate RNA/DNA cells prevent direct RNA
membership validation. Ploidy and DNA segmentation/clustering uncertainty remain.

## Executed counterexamples and primary reads

Both actually executed awk programs are extracted from frozen old/new launchers.
The fixture covers positions1,499999,500000,500001,999999,1000000,1000001, two
chromosomes/read groups and a missing-RG record. Explicit interval membership
provides independent expected answers. Each corrected program conserves all eight
eligible fragments and repairs three per-read assignments exactly. The net
positive bin-count difference is two because opposite changes in a shared bin
cancel; it must not be mistaken for the number of moved reads.

The primary trace queries the first20 chr1 500kb boundaries in each local archive:
SRR8117283, SRR8117281 and SRR8117282. All60 queries succeed. Independent filtering
finds22 eligible fragments:8 SNU-601,7 MKN-45 and7 NUGC-4. Ten start exactly at a
boundary:4,3 and3 respectively. Each old/new bin differs by one as predicted.
For example SNU-601 read496447 has SAM POS3000000, flag99 and MAPQ60; its correct
zero-based coordinate2999999 belongs to bin5, while the old expression gave6.
No genomic prevalence or total-error extrapolation is made from these loci.

Archive REFERENCE rows identify chr1 as NAME and SEQ_ID. The installed sam-dump
help declares one-based inclusive query endpoints. A known-read control connects
SRA zero-based REF_POS9991 with SAM POS9992 in both point and wider-window output.
Initial header-only loci therefore do not establish an assembly-name mismatch or
absent genome-wide coverage. SAM/SRA output, commands, errors, headers and hashes
are retained. Full archive checksums were not reassessed. The preserved original
SNU-601 BAM MD5 record agrees with its source table; that historical record is
not a new verification of the whole archive or FASTA identity.

All ten boundary reads carry sample-index read groups and lack CB/CR tags.
They cannot recover which original DNA cell received each read. The source BAM
and SRA representation retain different cell-tag information, so these reads
cannot establish a cell-level clone change. Independently regenerated counts
from suitable original material, followed by DNA-only resampling/segmentation
assessment, would establish whether the correction changes supported profiles.

## Verification, migration and remaining work

Probe23007104 completes in2m58s, process RSS148052KiB; reference metadata23007427
in20s; known-read control23007910 in59s. Independent audit23008067 completes in4s,
RSS83120KiB. It checks all60 query hashes, independently reapplies the read filters,
verifies read uniqueness and all22 fragment counts, and reconstructs the ten
misassignments by interval membership. Final fixtures23008603 complete in4s,
RSS79104KiB after the output-directory option is added. The final awk mapping
and all exact fixture requirements pass; shell syntax and diff checks pass.

The first probe23006938 failed its aggregate-count assertion before archive
queries. It incorrectly equated net bin changes with per-read moves. The repaired
check retains exact counts and verifies all per-read changes independently.
The failure, old/new outputs and every source version are preserved; no biological
criterion or numerical tolerance was relaxed.

Future corrected fetches can supply a fresh optional output directory:
`fetch_bins.sbatch LINE RUN OUT_BINS_DIR` or
`fetch_original_bam.sbatch LINE URL MD5 OUT_BINS_DIR`. Existing default paths and
overwrite guards remain. Pass any new count directory explicitly to downstream
DNA reconstruction and give resulting truth/scores new identities. Existing
count files do not become corrected merely because the source changes. No full
fetch, genome reprocessing, DNA re-clustering or RNA rescoring was run here.

The package wheel remains the previously validated351-Python/4-Rust artifact,
with one GPU-only skip. Package Python/Rust code is unchanged this round; the
full suite was not rerun for two shell fetchers. Evidence is under
`dev/andor_dna_coordinates_20261108/{probe_v1,probe_v2,reference_v1,reference_v2,
audit_v1,fixtures_final_v1}`. Historical DNA truth retains this additional
coordinate limitation. Original recovery, false-split/power and biological-null
gaps remain, with no SOTA or biological-validation promotion.

## Versioned primary-count repair submitted after the bounded audit

The verified upstream defect justifies regenerating the SNU-601 original BAM
counts: aggregated bins and the inspected SRA sample-index groups cannot recover
the affected original cell assignments. Job23009248 streams only that original
BAM into fresh `dev/andor_dna_coordinates_20261108/corrected_counts_v1/bins`,
requiring the published source MD5 before completing. It requests4 CPUs,16GB and
12h on innovation/io; live io partition maximum3d and QOS cpu500 limit were checked.
This is a submitted repair, not a completed or validated new DNA truth. Its
`rebuild_<job>/COMPLETE`, source MD5, output hashes, exit status and scheduler state
must be inspected before downstream reconstruction. No old directory is reused.

No improved clone result is presumed. Once counts are validated, reconstruct DNA
independently of RNA and assess coverage, DNA resampling and segmentation stability
before treating changed groups as stronger truth. Frozen evidence is unchanged;
the preservation record checks516 baseline files and69 unchanged Python/Rust/test
files for this round.

## Primary-count repair completed

Rebuild23009248 completes45m08s, source MD5 matches, and output hashes/COMPLETE
are recorded. Count audit23015231 passes exact conservation over405583339
fragments/1531 read groups, with832 autosomal boundary relocations in610 groups.
The original counts and all historical truth/scores remain unchanged. DNA
reconstruction/stability and independent numerical auditing are tracked in
[the follow-up protocol](ANDOR_DNA_STABILITY_2026-11-08.md); no clone improvement
or DNA-truth certification follows from count conservation alone.
