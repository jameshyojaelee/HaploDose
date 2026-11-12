# Synthetic RNA counting component before Bian measurement

This implements an exact gene-wise UMI counter behind existing RNA-native
alignment and gene-assignment tools. It tests synthetic reads only; no Bian
measurement or genomic outcome is opened. Passing can qualify one counting
component for a future frozen measurement protocol, not validate HaploDose
recovery,biological error rates or competitor performance.

## Prespecified numerical answer and pipeline

A fixed seed2026111201 produces a12000-base artificial contig. Five genes
include an overlapping pair and a two-exon gene with a canonical splice junction.
Fifteen60-base reads test same-UMI duplicates at different positions, one-base-
different UMIs, valid homopolymer UMIs, the same UMI in different genes, spliced
reads,intergenic/unmapped reads and a read mapping equally to two loci. The
hand-calculated molecule-count vector is G1=4,G2=2,G3=0,G4=0,G5=1. The exact
answer is saved before any tool runs. No numerical tolerance applies.

STAR2.7.11b maps the transcript reads. Existing Conda featureCounts assigns
unstranded exonic reads to genes and annotates assignments in BAM. Existing
Conda samtools exports the tagged SAM. The new counter requires a primary,
nonsupplementary uniquely mapped read with XS=Assigned,XN=1 and a single XT gene,
then counts exact unique UMIs within that gene. Same-UMI reads at different
positions collapse; the same UMI in distinct genes counts independently.
Missing assignment metadata,unknown barcode/name contracts or repeated primary
read IDs fail explicitly. All alignment records have an exclusion/count reason.

The [version-specific STAR parameters](https://github.com/alexdobin/STAR/blob/2.7.11b/source/parametersDefault)
describe basic STARsolo UMI filtering that removes homopolymers. That default is
not silently applied to the [source protocol's exact gene-wise UMI rule](https://link.springer.com/article/10.1186/s13059-018-1416-2).
The fixture deliberately retains valid homopolymer UMIs and distinguishes
one-base-different tags. This choice does not establish immunity to sequencing
errors or8-base UMI collisions. Biological molecule counts and their sampling
law still require qualification; technically exact deduplication alone does not
make an NB uncertainty claim valid.

A separate hand-constructed tagged-SAM fixture checks the counter independently
of alignment and featureCounts, including rejected evidence. The end-to-end
check must reproduce every positive read's intended gene and all five exact
gene counts. A failed alignment,assignment or counter check is preserved with
its failing exit status. There is no post-failure relaxation of expected counts.

## Compute and scientific scope

Job23231315:nslab account/nslab QOS,CPU2CPU4GB1h,no array. Existing STAR module
and project Conda tools,isolated Python runtime,neutral directory and141 source
identities before/after. AST and shell syntax pass; executable result pending
at submission. Exact commands,tool versions/hashes,read identities,expected
counts and rejection ledger are retained. No new dependency or shared install.

The real Bian adapter trimming,barcode/UMI handling,gene counts,DNA fragments,
callability,normal references and genomic interpretation are not validated by
this fixture. The original candidate/reference/comparator/scorer measurement
freeze still applies. Current recovery,false-split/power,uncertainty and all-
competitor criteria remain unchanged. Biological n=0. No production or frozen
result changes.

[Frozen design](../dev/bian_count_fixture_20261112/DESIGN.md),
[exact gene/UMI counter](../dev/bian_count_fixture_20261112/count.py),
[synthetic pipeline and independent expected result](../dev/bian_count_fixture_20261112/run.py),
[source identities](../dev/bian_count_fixture_20261112/SNAPSHOT.json),
[submission](../dev/bian_count_fixture_20261112/submission.json),
[qualified shared-reference inputs](BIAN_REFERENCE_CONTRACT_2026-11-12.md).

## Completed counting result

Job23231315 completed0:0 on ne1dc6-003 in4s. The timed process took1.65s with
peak260676KiB; Slurm did not sample useful RSS. All141 source identities pass
before and after. The manual counter fixtures and complete synthetic alignment/
gene-assignment/deduplication chain pass without changes to the expected counts.

The output is exactly G1=4,G2=2,G3=0,G4=0,G5=1:seven gene-specific distinct UMIs
from11 assigned reads. All15 input reads retain unique primary identities; one
extra secondary alignment is recorded separately. The four excluded primary
reads are ambiguous between genes,intergenic,nonuniquely mapped and unmapped.
The spliced gene's two reads count once. Homopolymer UMIs and one-base-different
UMIs remain distinct; the same UMI in two genes counts in each gene. Every
expected positive read maps to its intended gene.

This is a passed computational measurement-component check. It does not validate
real Bian tag extraction,adapter removal,reference matching,DNA sensitivity,
RNA/DNA dosage or clone identity. It does not repair the outstanding recovery,
state-null/power or uncertainty failures or change any package defaults.
The future full measurement protocol must preserve these count semantics and
specify its real-read preprocessing before it can open genomic outcomes.

[Exact expected/actual result](../dev/bian_count_fixture_20261112/run_23231315/results/REPORT.json),
[full read exclusion/count ledger](../dev/bian_count_fixture_20261112/run_23231315/results/COUNTS.json),
[independent manual fixtures](../dev/bian_count_fixture_20261112/run_23231315/results/MANUAL_FIXTURE.json),
[executed tool identities](../dev/bian_count_fixture_20261112/run_23231315/results/TOOLS.json).
