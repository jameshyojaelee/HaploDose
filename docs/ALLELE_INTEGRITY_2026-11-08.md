# One locked Kinker line retains donor-specific counts and phase

The lexical diagnostic selects `2313287_STOMACH` from189 artifact-bearing lines
before inspecting allele outcomes. The full enumeration and selection rule are
saved in `dev/kinker_allele_integrity_20261108/selection.json`. This is one
historical mechanism investigation, not an estimate of population performance.

Compute job23000714 completes in41s with312992KiB peak RSS. Pool7_16 contains54
selected cells and11415 retained SNPs; Pool8_16 contains50 cells and10882 retained
SNPs. No selected barcodes are missing. Source barcode maps are unique and donor
specific. Independently selecting these columns from the original channel AD/DP
matrices reproduces the saved donor matrices exactly. Reconstructed donor totals
also reproduce the VCF AD/DP/OTH fields and the retained SNP identities/order.

Across22 chromosomes,3500 genotype/phase records retain the selected donor
identity, variant identities and unordered genotype. All24576 exported cell-SNP
rows agree with the corresponding donor pileup counts and phased heterozygous
orientation. These equalities use exact integer counts and string identities,
without a numerical tolerance or a new filtering rule. No allele artifacts or
historical comparator outputs were rewritten, and no inference was rerun.

This does not establish the original reference assembly or the exact historical
source version. The current phase script is captured as context; matching its
predictions is not proof that this exact source was executed historically. The
check validates correspondence of retained rows, not completeness of every
scientifically eligible SNP, phase accuracy against DNA, absence of mapping bias,
or biological LOH. No read inspection is triggered because this bounded count
and orientation trace finds no discrepancy. Independent DNA alleles and reserved
SNP evidence remain necessary for stronger biological claims.

Implementation: `benchmarks/hardening/kinker_allele_integrity.py`. Immutable
execution inputs, environment, output hashes and report are in
`dev/kinker_allele_integrity_20261108/{source_v1,audit_v1,analysis_23000714}`.

## Bounded original-BAM header follow-up

Job23000741 completes in3s using SAMtools1.19.2. Headers from SRR12557907
(Pool7_16) and SRR12557908 (Pool8_16) both record STAR2.5.1b alignment against
`refdata-cellranger-hg19-1.2.0/star`; chromosome1 has length249250621 in each.
This supports consistency with the named hg19 phasing panel and narrows the
assembly gap reported above. Reference FASTA checksum identity and the exact
historical phasing execution remain unestablished. The header command adds its
own samtools PG entry; that new entry is inspection provenance, not an original
alignment command. No reads are reprocessed or inspected at loci, and no full
BAM checksum is claimed. Headers, hashes, source, job record and file metadata
are retained under `dev/kinker_allele_integrity_20261108/headers_v1`.
