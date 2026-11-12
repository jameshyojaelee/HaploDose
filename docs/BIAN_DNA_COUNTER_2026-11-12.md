# Independent DNA fragment counter before truth construction

Implemented an isolated strict paired-DNA counter for the future native-hg19
measurement path. It rejects nonreciprocal mate coordinates and falsely declared
query-name ordering,checks CIGAR query/reference spans and retains complete empty/
rejection accounting. Synthetic fixtures pass;no specimen DNA was processed.
This is computational preparation,not CNA truth or a performance improvement.

The estimand is one fragment-midpoint observation per accepted template under a
declared duplicate-removal rule. The counter does not prove deduplication or
independent molecules. Proper-pair flags and mapping quality are measurements;
the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) defines
mate-coordinate,reference-span and Phred conventions. The nominal MAPQ20 screen
corresponds to1% error only under calibrated mapper qualities;255 is unavailable,
not the highest confidence. This does not establish specimen false-split control.

The counter validates both primary mates,reciprocal RNEXT/PNEXT and orientation
metadata,CIGAR lengths,reference bounds and the exact SAM reference dictionary.
It rejects N reference skips for this genomic Bismark/Bowtie2 path,whose contract
does not include spliced RNA alignments. It verifies actual lexicographic query
ordering with current-group storage rather than trusting SO:queryname. Autosome
counts and all other contigs remain separate. More than two primary records are
rejected with bounded storage. Empty outputs retain an explicit no-evidence flag.

Job23232108:nslab CPU1CPU1GB1h,COMPLETED0:0 in1s;timed process0.54s,peak11772KiB.
All174 identities pass before and after. The independent15-group/32-alignment
fixture has4 accepted autosomal templates,one other-contig template and10 rejected
groups. Exact Fraction midpoint calculations give chr1 bins0/1 counts3/1 and
chrX bin0 count1. Low/unknown MAPQ on either mate,flags,extra primary alignments,
reciprocal-coordinate defects,spliced DNA,malformed CIGAR,query-length/bounds,
false sorting and dictionary mismatch behave as prespecified. Empty input is
accounted for. The native reference assets were qualified separately.

Only pure constants/functions from the existing CRC classifier were executed
through an isolated AST,without importing its I/O helpers. It accepts the three
constructed bad RNEXT,PNEXT and spliced-DNA cases;the new classifier rejects them.
This establishes an input-contract gap,not occurrence of those faults in past
biological data. No CRC source,shared installation,historical result or public
HaploDose code was changed,and no previous scores were recalculated.

Independent genomic validation still needs qualified aligned/deduplicated DNA,
callable support,genomic-gap handling,normal/reference variation,resampling and
resolution-matched comparisons. Bin width is a counting convention,not validated
assay resolution. No purity,ploidy,LOH,clone or uncertainty claim follows. The
full measurement/evaluation freeze,secondary/descriptive role and original
recovery,false-split/power,uncertainty and competitor criteria remain unchanged.

[Frozen design](../dev/bian_dna_counter_20261112/DESIGN.md),
[strict counter](../dev/bian_dna_counter_20261112/fragments.py),
[independent fixture and pure legacy comparison](../dev/bian_dna_counter_20261112/run.py),
[exact result](../dev/bian_dna_counter_20261112/run_23232108/results/REPORT.json),
[source identities](../dev/bian_dna_counter_20261112/SNAPSHOT.json).
