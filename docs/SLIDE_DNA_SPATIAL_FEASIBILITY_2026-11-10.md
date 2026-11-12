# Colorectal serial-section validation feasibility

The Figure 4 DNA and RNA libraries now have verified raw-run identities,
listed bead-level count files, and a reconstructable published coordinate
map. The pair can support a regional validation pilot once registration
provenance, DNA coverage, and reference choice are qualified. This step
evaluated geometry only. It produced no RNA-versus-DNA accuracy estimate.

## Exact measurement chain

All five section aliases match both ENA sample and experiment aliases
uniquely among the 18 runs of PRJNA768453. Figure 4 uses adjacent sections
of one colorectal tumour. The Figure 3 slide samples a different region,
and its single-cell DNA libraries sample the wider tumour.

| Role | Section | ENA sample | Run |
|---|---|---|---|
| Figure 4 spatial DNA | `human_colon_cancer_4_dna_200114_13` | SAMN22045091 | SRR16203709 |
| Figure 4 spatial RNA | `human_colon_cancer_4_rna_200102_06` | SAMN22045092 | SRR16203708 |
| Figure 3 spatial DNA | `human_colon_cancer_3_dna_191204_19` | SAMN22045088 | SRR16203712 |
| Figure 3 single-cell DNA | `human_colon_cancer_3_scwgs_200227_01` | SAMN22045089 | SRR16203711 |
| Figure 3 single-cell DNA | `human_colon_cancer_3_scwgs_210402_02` | SAMN22045090 | SRR16203710 |

The portal lists 28 files, including the Figure 4 DNA sparse 1 Mb counts
(117,216,767 bytes), RNA sparse expression (162,991,889 bytes), both bead
coordinate files, `hg19_1Mb_bins.txt`, and `human_genes.txt`. The DNA source
also reads `sparse_frags_1Mb.txt`, which is absent from that inventory.
Raw paired FASTQs total 13,461,610,657 bytes for Figure 4 DNA and
9,967,660,597 bytes for Figure 4 RNA. A separate normal-colon WGS run,
SRR16203706, is listed. Its donor correspondence still requires evidence.
Portal authentication may be needed to retrieve the listed counts. The
Zenodo record lists a 26,335,256-byte archive and an MD5. Only its record
was read. The author repository lists the registration file and separate
clone labels. Pinned RNA and DNA analysis scripts matched their Git blob
hashes and were read as source.

[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9301586/),
[section manifest](https://singlecell.broadinstitute.org/single_cell/study/SCP1278/spatial-genomics-enables-multi-modal-study-of-clonal-heterogeneity-in-tissues),
[repository at the verified commit](https://github.com/buenrostrolab/slide_dna_seq_analysis/tree/41142f767075df4ae84d2142cf4c771317a5da2f),
[archive record](https://doi.org/10.5281/zenodo.5553305).

## Coordinate reconstruction

The RNA script scales its two coordinate columns by 0.65, keeps points
less than 1,600 micrometres from their mean, and reads the registered
coordinate file in that row order. There are 34,263 unique RNA barcodes,
34,233 retained positions, exactly 34,233 registered rows with two
columns, and 40,766 unique DNA barcodes. All three files matched the Git
tree hashes at `41142f767075df4ae84d2142cf4c771317a5da2f`.

Fit a 2D affine map on even-numbered array indices, then measure Euclidean
errors on the held-out odd indices. The fitted linear part has determinant
0.99999996. Median and 99th-percentile errors are 0.03887 and 0.06529 in
registered coordinate units. Reversing the registration rows raises them
to 1,049.94 and 1,493.46. This reproduces the coordinate mapping and tests
row order independently of expression, DNA coverage, or clone labels.
The fit residual measures coordinate export agreement. Tissue alignment
error and the biological landmarks used to choose the transform remain
unverified. The published RNA filter explicitly uses the DNA-derived
subclone mask and is therefore ineligible as RNA-only input selection.

## What the pair could establish

Freeze a region definition using morphology or independently documented
geometry, and record landmarks or held-out image overlap and displacement
error. Keep RNA method inputs and predictions separate from DNA-derived
masks and clone labels. Match hg19 genomic support before scoring. Build
DNA-only regional relative coverage from bead counts with a justified
baseline and a callable mask. Use independent fragment partitions where
available, and test sensitivity to bead aggregation, thinning and
displacements supported by the geometry. Compare RNA profile predictions
with matched pooled RNA ratios and a shared smoothing rule on the same
regions and bins. Recheck conclusions across credible DNA baselines and
alignments. Record uncovered regions and registration-dependent outcomes.
Start with existing bead counts. Only an unresolved count, coverage, or
fragment-provenance discrepancy justifies inspecting bounded reads.

One tumour permits a regional compatibility case study. Prospective
population superiority and individual spot or cell membership require
additional specimens and measurement correspondence. The registered
reference candidate also requires an independent RNA reference field with
plausible state transport, which the present sources have yet to establish.
RNA and DNA count files and clone labels remained unopened.

## Compute and provenance

All jobs requested one CPU, 4 GB and 30 minutes on CPU node003, with an
isolated wheel and immutable source snapshots. Each completed with exit
0:0. The status of each registry request is recorded separately from
scheduler completion. The first follow-up received HTTP 406 from the
portal. The explicit JSON `Accept` header recovered its inventory. The
PMC supplement requests returned non-PDF responses, the Springer CDN
hostname did not resolve, and the Europe PMC full-text requests returned
404 and 500. Registration protocol verification remains open.

| Job | Account | Work | Elapsed | Process peak KiB |
|---|---|---|---:|---:|
| 23186010 | innovation | ENA, Git tree, Zenodo | 7 s | 20,752 |
| 23186036 | nslab | First portal and protocol fetch, pinned scripts | 5 s | 20,876 |
| 23186047 | nslab | Portal inventory recovery and protocol retry | 5 s | 20,932 |
| 23186060 | innovation | Hash-verified coordinate check | 7 s | 71,416 |

[Design](../dev/spatial_slide_dna_metadata_20261110/DESIGN.md),
[registry report](../dev/spatial_slide_dna_metadata_20261110/metadata_v1/report.json),
[first follow-up](../dev/spatial_slide_dna_metadata_20261110/followup_v1/report.json),
[retry report](../dev/spatial_slide_dna_metadata_20261110/followup_v2/report.json),
[public file inventory](../dev/spatial_slide_dna_metadata_20261110/followup_v2/scp1278.json),
[geometry result](../dev/spatial_slide_dna_metadata_20261110/geometry_v1/report.json),
[Slurm accounting](../dev/spatial_slide_dna_metadata_20261110/accounting_full.tsv).
