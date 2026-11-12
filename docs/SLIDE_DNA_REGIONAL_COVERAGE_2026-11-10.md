# Regional DNA coverage in the colorectal serial section

Four fixed archive intervals yielded 1,432,964 coordinate-deduplicated,
autosomal DNA fragments inside the bead puck. The two disjoint fragment
halves support most 10 and 20 Mb tile-versus-rest comparisons under the
prespecified count planning rule. Support at 1 Mb is sparse. The first half
has zero windows exceeding the prespecified large-effect threshold among
windows with sufficient counts, so the sign-replication test has zero
eligible effects. This one tumour supplies technical coverage evidence.
Physical section registration and RNA-to-DNA accuracy remain unresolved.

## Preserved measurement chain

The first 2,048 archive spots replayed the published bead correspondence
for 815 pairs. Native SRA Phred qualities yielded 214 canonical aligned
rows, 207 unique, compared with 211 and 204 using the earlier ENA export.
One fragment row appeared only in the original alignment and four only
with native qualities. The per-base comparison found all 28,525 qualities
per mate differed. ENA assigned a constant Phred 30, except for five mate
2 reads assigned Phred 3. Applying the documented SRA Lite rejection rule
to the native prefix reproduced both original uncompressed FASTQ SHA-256
digests and all 211 aligned rows by spot ID. The saved native and original
files differ in alignment output order because the replay used `--reorder`.
The main extraction uses native per-base qualities.
[Native comparison](../dev/spatial_slide_dna_coverage_20261110/control_native_report.json),
[base-quality comparison](../dev/spatial_slide_dna_coverage_20261110/quality_23186255/report.json),
[identity-ordered Lite replay](../dev/spatial_slide_dna_coverage_20261110/control_lite_identity_report.json),
[NCBI quality format](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/).

The four 4,194,304-spot intervals start at 1, 80,000,001, 160,000,001 and
240,000,001 in SRR16203709. Across 16,777,216 spots, 6,625,299 barcodes
match a published bead after reverse complementation. Bowtie2 alignment
against the checksum-verified hg19 reference and proper-pair MAPQ ≥30
filtering yield 1,914,356 paired groups, including 1,874,059 canonical
fragment rows. Coordinate deduplication within each bead leaves 1,481,406
unique fragments. The in-puck autosomal restriction leaves 1,432,964
fragments across 40,331 beads. The seeded technical split yields 716,032
and 716,932 disjoint fragments. Each extraction records FASTQ and fragment
hashes, exact commands, modules and resource logs.
[Locked design](../dev/spatial_slide_dna_coverage_20261110/DESIGN.md),
[extraction commands](../dev/spatial_slide_dna_coverage_20261110/jobs_extraction.jsonl),
[input identity and full regional report](../dev/spatial_slide_dna_coverage_20261110/coverage_23187918/results/report.json).

## Resolution and boundary sensitivity

Nine fixed 1 mm tiles are compared with the rest of the puck. The endpoint
is the tile-versus-rest log2 DNA fragment-proportion ratio over eligible
autosomal windows. Windows crossing major hg19 assembly gaps are masked.
The count planning threshold is a Poisson delta half-width of 0.292481
log2 units, half the pure 3:2 effect. It corresponds to 187 fragments in
each side of a balanced window in each half. This planning calculation
assumes Poisson sampling. Biological variability and regional technical
bias remain outside that approximation. The same-assay rest-of-puck
baseline has uncertain tumour composition and ploidy.

| Width | Eligible tile-windows | Both halves meet count rule | Half disagreement ≤0.292481 among supported | Largest tile p95 half disagreement | Largest tile p95 change after 50 µm shift |
|---:|---:|---:|---:|---:|---:|
| 1 Mb | 23,958 | 32 | 30/32 | 0.443 | 0.116 |
| 5 Mb | 4,464 | 3,133 | 2,897/3,133 | 0.368 | 0.099 |
| 10 Mb | 2,043 | 1,867 | 1,812/1,867 | 0.307 | 0.084 |
| 20 Mb | 819 | 810 | 803/810 | 0.294 | 0.081 |

The weakest tile has 5/496 supported 5 Mb windows and 61/227 supported
10 Mb windows. Every tile has 90/91 supported 20 Mb windows. One fixed
within-puck bead permutation produces 13, 3,172, 2,031 and 810 windows
meeting the same count rule at the four respective resolutions. Count
support and technical agreement therefore cannot establish spatial CNA
contrast by themselves. Both the registered tiles and permutation have
zero first-half supported windows with |log2 ratio| ≥ log2(1.5). The
50 µm shifts are specified perturbations. The physical registration error
between sections still requires independent landmarks. A focal or mixed
alteration can be diluted by tile averaging or shared with the rest of
the puck. These results do not support scoring spot genotypes, or choosing
new tiles based on an RNA method's output. A future regional RNA–DNA
comparison needs independently chosen tissue regions, a registration
error budget, and demonstrated DNA sensitivity to the intended contrasts.
No RNA counts or clone labels were opened in this analysis.
[Complete windows including gaps and failed support](../dev/spatial_slide_dna_coverage_20261110/coverage_23187918/results/profiles.tsv.gz).

## Compute and failure recovery

| Job | Account | CPUs / RAM / time limit | Result |
|---|---|---|---|
| 23186240 | innovation | 2 / 12 GB / 30 min | Native quality replay, 0:0, 6 s |
| 23186255 | nslab | 1 / 4 GB / 30 min | Base-quality probe, 0:0, 4 s |
| 23186267 | nslab | 2 / 12 GB / 30 min | SRA Lite replay, 0:0, 6 s |
| 23186817, 23186818, 23186819, 23186820 | innovation, nslab, innovation, nslab | 4 / 12 GB / 2 h each | All 0:0, 2m14s–2m24s, max batch RSS 3,355,016 KiB |
| 23187150 | nslab | 1 / 8 GB / 1 h | 1:0, stopped on the checksum working directory before analysis |
| 23187918 | nslab | 1 / 8 GB / 1 h | 0:0, 15 s, process peak 185,208 KiB |

The scorer retry checked the relative design paths from the repository
directory before entering a fresh output directory. It verified the same
design, source and original launcher hashes. The numerical fixtures verify
analytic log ratios, Poisson planning widths, gap boundaries and fragment
conservation. The running extraction jobs were left alone when the account
instruction changed. The retry and later submissions use nslab.
[Failed launch log](../dev/spatial_slide_dna_coverage_20261110/coverage_23187150.log),
[retry command](../dev/spatial_slide_dna_coverage_20261110/jobs_score_v2.json),
[input hashes](../dev/spatial_slide_dna_coverage_20261110/coverage_23187918/inputs.sha256),
[process resources](../dev/spatial_slide_dna_coverage_20261110/coverage_23187918/resources.txt),
[Slurm accounting](../dev/spatial_slide_dna_coverage_20261110/accounting_complete.tsv).
