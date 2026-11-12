# Prospective truth construction: source review

This review narrows the outstanding measurement requirements for the two existing
public leads. It does not accept a cohort or supply biological replication. Only
protocol text and source code were inspected; no per-cell measurements, clone
tables, fragments or raw reads were opened. Author-reported findings appearing
alongside the methods were visible, so this is not a fully blinded literature review.
Current repository and inspected blob identities are recorded in
`dev/prospective_source_review_20261108/metadata.json`. No external code was run.
Current public code is not proof of the version executed for published artifacts.

## Hybrid BAG-seq

The protocol maps to hg19, uses the same BAG identity for DNA/RNA templates and
separates molecular layers using annotation overlap. Its final DNA filter excludes
RNA hotspots derived from RNA-only libraries of two normal and five tumour
tissues. DNA-only controls are available. Published copy-number clustering uses
300 empirical bins, spanning 4.7–36.5 Mb, with a median 9.5 Mb. Adjacent-bin
summaries enter clustering; segmentation is for aggregate profile displays.
[Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12476631/)

The inspected [DNA clustering script](https://github.com/siranli01/DNA_RNA/blob/c67e399041a453d048b345b129d0f2534b05a830/Clusteing.based.on.copy.number.with.sliding.window.R)
accepts a normalized bin matrix and genomic boundaries. It augments that matrix
with two/three-bin medians, removes windows crossing chromosomes, then uses PCA,
neighbours and Seurat clustering. It does not take RNA groups as an explicit
input. This establishes the visible script's input structure, not full independence
of the preprocessing, supplied matrices or author decisions. The upstream
[pipeline README](https://github.com/levycshl/bag-pipe/blob/e9acdf729afd0285b5a9a2d9608227b74198a142/README.md)
describes barcode extraction, HISAT2 mapping and template summaries; an exact
archived run configuration and the final layer/mask construction are still needed.

[Inference] The admission test must distinguish DNA-specific CNA evidence from
RNA contamination or filtering-induced agreement. Freeze an independently
specified mask, reconstruct DNA groups without RNA group labels, and check their
stability against the corresponding DNA-only protocol. If the published mask is
used, report its RNA-derived provenance and test a separately derived control;
do not describe the resulting measurements as completely independent of RNA.
Evaluate only contrasts resolved by DNA at its empirical bins and coverage.
The published bin size does not support automatically calling 1-Mb truth reliable.

The smallest decisive check is a predeclared matching table linking each hybrid
library, donor, DNA-only control, RNA-only control, barcode convention, assembly,
bin table and mask version. Missing joins remain failures of correspondence.
After a frozen evaluation design, DNA resampling and coverage reduction must
establish which contrasts and memberships are actually resolvable. Agreement of
RNA with published labels alone would not pass this gate.

## DEFND-seq

The methods describe initial cell selection from RNA UMI ranks, whitelist and
RNA/DNA barcode concordance, and filtering by both RNA and DNA counts. Human
RNA processing uses GRCh38/Gencode v24 outside the species-mixture experiment.
The assay separates RNA and DNA libraries but acknowledges residual accessibility
bias. These are assay and selection properties, not an artifact-specific barcode
mapping or a validated DNA sensitivity estimate.
[Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12272635/)

The inspected [DNA driver](https://github.com/simslab/dna10x/blob/2948ac7c31e08a3d46b2d5c5aa9cd34507d038d9/dna10x.py)
takes a reference, barcode whitelist, barcode position, insert-size/alignment-score
settings and an optional whitelist reverse-complement flag. The
[barcode correction](https://github.com/simslab/dna10x/blob/2948ac7c31e08a3d46b2d5c5aa9cd34507d038d9/error_correct.py)
uses sequence quality and whitelist barcode frequencies. None of these settings
should be guessed from a GEO library name.

The [fragment writer](https://github.com/simslab/dna10x/blob/2948ac7c31e08a3d46b2d5c5aa9cd34507d038d9/count_dna.py)
groups records by genomic endpoints and corrected barcode, writes one row per
such key and places its multiplicity in column five. It leaves reference-hash
header fields empty. [Inference] If this version produced the deposited files,
summing column five as independent DNA fragments would overcount duplicated
records. That provenance must be established before selecting a counting rule;
neither a pipeline-name header nor code availability establishes assembly identity.

The smallest admissible check after evaluation freezing is a bounded fragment
header/coordinate and barcode crosswalk inspection tied to the published RNA
matrix and archived command settings. Require one-to-one library-qualified
identity, declared coordinate conventions and explicit handling of duplicate
fragment multiplicities. If provenance disagrees, use bounded raw alignment/read
inspection where available. An unavailable raw DNA chain must be reported as
unverified, not repaired by assuming the RNA assembly or barcode spelling.

## Decision

Both studies remain possible same-cell challenges, not accepted validation panels.
This review supplies concrete gates rather than clearing them. No threshold,
cohort, specimen or resolution was selected using HaploDose outcomes. There is
still no verified panel of 59 independent single-CNA null specimens; cells,
libraries, mixed-donor pools and repeated computational draws cannot replace it.
Former evaluation cohorts remain regression-only and Tahoe remains withdrawn.
