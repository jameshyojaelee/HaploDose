# Recover colorectal spatial DNA bead identity from the raw archive

The processed matrices require Single Cell Portal sign-in through its
download API. The normalized SRA object for Figure 4 DNA preserves the
14-base technical read needed to recover spatial bead identities. In a
bounded prefix, reverse complementation links 815 of 2,048 sequenced
spots to the published coordinates. Bowtie2 and the independent filter
check yield 204 coordinate-unique fragments. This establishes a route
to independent regional DNA counts using the public primary archive.

## Access and primary measurement chain

The public SCP1278 metadata inventory is available. Job 23186098 tests
the two exact count files through the API and legacy download routes.
HEAD and ranged GET return HTTP 401 at the API and HTTP 500 at the
legacy route. The probe does not follow redirects or read response
bodies. The NLM public cloud lists the author manuscript PMC9301586.1
with metadata, XML and text only. The protocol supplement is absent from
that inventory. GitHub API rate limiting interrupts the independent
benchmark source lookup. The later `git ls-remote` and raw source URLs
recover the benchmark documentation at a pinned commit.

The published DNA alignment template starts from three original read
streams, extracts the barcode from the middle stream, aligns the outer
paired reads with Bowtie2 and filters proper pairs at MAPQ 30. Its
sparse-count script reverse-complements the barcode from the coordinate
table. The ENA record for SRR16203709 and SRX12488063 lists two public
FASTQs and an empty submitted-file field. Job 23186120 downloads exactly
the first 1 MiB compressed per mate and checks 2,048 synchronized pairs.
All sampled mates are 35 bases. Neither mate contains the expected
linker in either orientation. Exact 14-base coordinate matches at
arbitrary positions occur in 6–11 reads per mate/orientation. These
matches alone provide no evidence for a barcode layout.

Job 23186130 reads only `READ_START`, `READ_LEN`, and `READ_TYPE` from the
normalized SRA object via `vdb-dump` 3.2.1. All 128 rows give starts
0/35/49 and lengths 35/14/35, with the middle read typed as technical.
Job 23186147 extracts the 14-base middle string on 2,048 spots, checks
against all 40,766 exact coordinate barcodes and stores aggregate counts.
The reverse complement matches 815 spots. Direct orientation and a
one-base rotation of the reverse complement match zero. The two genomic
mate prefixes yield zero and two matches. There were no ambiguous bases
in the middle reads. These are prefix diagnostics and do not estimate
whole-run recovery, DNA callability, or biological replication.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9301586/),
[public study inventory](https://singlecell.broadinstitute.org/single_cell/study/SCP1278),
[SRA run](https://www.ncbi.nlm.nih.gov/sra/SRR16203709),
[authors' alignment template](https://github.com/buenrostrolab/slide_dna_seq_analysis/blob/41142f767075df4ae84d2142cf4c771317a5da2f/scripts/run/slidedna_alignment_template.sh),
[authors' count writer](https://github.com/buenrostrolab/slide_dna_seq_analysis/blob/41142f767075df4ae84d2142cf4c771317a5da2f/scripts/preprocess/get_sparse_counts.py),
[NLM cloud inventory documentation](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/).

## Spatial comparators

The ST-CNABench source at `e0112d1a67572c95ecc1b227c72de7fef575c070`
confirms that its colorectal tutorial uses PRJNA768453 and expects an
already prepared spatially binned `SlideRNAseq_human_colon_bin10` input
with barcode mappings, clone labels and clone-level CNA reference
files. The tutorial identifies the public raw source and gives no
processed colorectal bundle location. Those historical scores would
need independently verified inputs and a matched resolution before
they could support a head-to-head comparison.

The same source gives a separate public cSCC demo for `P6_vis_rep1`
and names matched tumour/normal WES coverage and FACETS files.
GEO identifies this Visium section as GSM4565825, processed with
Space Ranger against GRCh38. This is a candidate case for aggregate
DNA comparison. One patient's WES has no spot-level genotype labels.
The archive and scores were left unopened.

[ST-CNABench paper](https://www.nature.com/articles/s41467-026-77500-5),
[colorectal comparator tutorial](https://github.com/YangLabHKUST/ST-CNABench/blob/e0112d1a67572c95ecc1b227c72de7fef575c070/docs/tutorials/subclone_identification_slidednaseq.md),
[cSCC configuration](https://github.com/YangLabHKUST/ST-CNABench/blob/e0112d1a67572c95ecc1b227c72de7fef575c070/configs/examples/cscc_demo/data.yaml),
[cSCC Visium record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4565825).

## Compute record

| Job | Account and partition | Requested CPU / RAM / wall | Work |
|---|---|---|---|
| 23186098 | nslab CPU | 1 / 4 GB / 30 min | Access and publication inventory |
| 23186120 | innovation CPU | 1 / 4 GB / 30 min | ENA mate prefix and comparator source |
| 23186130 | nslab CPU | 1 / 4 GB / 30 min | Normalized SRA layout |
| 23186147 | nslab CPU | 1 / 4 GB / 30 min | Exact technical barcode match |
| 23186180 | innovation I/O | 1 / 8 GB / 1 h | Checksum-verified hg19 reference staging |
| 23186193 | nslab CPU | 2 / 12 GB / 1 h | Bounded cross-archive alignment |

The first four jobs finish with exit 0:0 in 6, 38, 4 and 4 seconds on
CPU node003. Process peaks are 21,304, 30,980, 74,092 and 77,504 KiB.
Index staging finishes 0:0 in 2m31s on I/O node002 using the UCSC hg19
p13 no-alt analysis set. The publisher's MD5 matches and the extracted
reference hashes and chromosome dictionary are recorded. The authors'
exact index build remains unknown. The alignment job finishes 0:0 in
13s on CPU node003, with Bowtie2 peak 3,295,284 KiB. Both VDB genomic
sequences match the ENA sequences and spot IDs exactly for all 2,048
records. ENA qualities and the VDB technical barcode were used for the
815 bead-linked pairs. SAMtools and an independent Python filter agree
on 213 complete proper pairs with MAPQ at least 30 in each mate. Of
these, 211 are on canonical chromosomes and 204 remain after exact
bead/chromosome/start/template-length deduplication. They cover 201
beads and 197 nominal 1 Mb bins. These counts qualify the recovery
path. They cannot estimate callable regions or substitute for the
authors' full barcode and duplicate handling. Before scoring regional
accuracy, obtain a defensible registration uncertainty estimate and
independently lock the RNA and DNA inputs and callable support.

[Access report](../dev/spatial_slide_dna_access_20261110/access_v1/report.json),
[bounded mate and benchmark source report](../dev/spatial_slide_dna_access_20261110/reads_v1/report.json),
[SRA layout](../dev/spatial_slide_dna_access_20261110/sra_layout_23186130/layout.tsv),
[technical barcode report](../dev/spatial_slide_dna_access_20261110/sra_barcode_23186147/report.json),
[bounded alignment counts](../dev/spatial_slide_dna_access_20261110/alignment_23186193/alignment.json),
[Slurm accounting](../dev/spatial_slide_dna_access_20261110/accounting_complete.tsv),
[UCSC reference checksums](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/md5sum.txt).
