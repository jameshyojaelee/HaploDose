# Bulk output preserves direct gene RNA effects

Bulk runs now retain direct gene effects before intermediate bin aggregation,
alongside the existing binned ratios and HMM output. This implements the supported
information-preservation step from the [matched gene study](GENE_AGGREGATION_2026-11-08.md).
It does not introduce a superior DNA estimator: the stronger naive comparator
remains unbeaten, and no genomic accuracy or uncertainty claim is promoted.

## Scientific and measurement contract

The output describes per-gene normalized RNA abundance relative to pooled normal
raw gene rates, conditional on the retained genes and declared reference samples.
For gene g and sample s, the statistic is
`log((G_gs+.5)/(L_s*mu_g+.5))-centre_s`, where L sums retained input gene counts,
mu is the pooled normal gene rate, and centre is the median of within-bin
gene-effect medians. The half-count offset is preserved from the established
development path rather than selected on DNA agreement. It is a normalization
choice, not a biological detection threshold.

Each normal excludes itself from its reference counts and library. Feature
selection still uses the common all-normal universe; this is explicitly not
complete-workflow leave-one-out calibration. The data may be reads or UMIs
according to the declared raw input. Programme effects, HMM dosage states,
smoothing and conditional SEs are not transferred to gene effects.

The strongest alternative to a DNA-driven gene effect remains transcription,
reference composition or count-scale mismatch. Gene retention prevents an
aggregation loss; it does not distinguish those causes. Corresponding DNA at
matched genomic resolution and independent RNA replication would establish
whether a preserved regional signal reflects dosage. Existing exposed cohort
comparisons remain development evidence.

## Interface and resource behavior

`gene_rna_log_ratio.npy` is float64 with rows in `samples.tsv` order and columns
in `genes.tsv.gz` order. The latter retains genomic coordinates and original
input-column correspondence. `gene_ratio_metadata.tsv` records each sample's
mapped input count, reference input count, reference sample number and centering.
`run.json` records the formula, axes, precision, reference/feature-selection
contract, units and unvalidated interpretation. The run manifest hashes both
new artifacts. Historical readers and frozen files remain unchanged.

Output is streamed through a memory map; per-bin gene memberships are computed
once. Additional disk storage is approximately eight bytes per retained
gene/sample pair, plus the NPY header and metadata. No additional full dense
gene/sample matrix is allocated. The existing production HMM and binned profile
calculations remain unchanged.

## Validation gates

The focal counterexample has two three-gene bins whose median effects are zero,
but contains retained gene effects log15 and log31. It demonstrates information
lost by nested medians without introducing a convenient tolerance. Public-run
tests independently sum scalar reference counts under deliberately shuffled
gene/sample order, including literal sample ID `001` and each normal's exclusion.
Both sum and median modes must produce the same direct raw gene effects.

The independent parity analysis executes the preserved old wheel and new wheel
on the same raw fixture in both modes. All six old scientific artifacts must
match exactly, comparing decompressed gzip payloads; old run metadata must match
after removing only the new gene-export declaration. The scalar calculation
must match every emitted gene effect exactly. A separate replay requires package
arithmetic to reproduce the frozen direct-gene reference in all nine development
refits. These checks assess computation and identity, not independent superiority.

## Executed evidence and decision

Environment recreation 23004347 completes in 3m19s, using existing dependencies
in a new isolated prefix. Build/test 23004495 completes in 2m25s: 351 Python
tests and four Rust tests pass, with one GPU-only skip. The full test process
peaks at 2527524 KiB RSS. All preexisting tests remain unchanged.

Parity job 23005012 completes in 26s. Both bulk modes preserve all six original
scientific artifacts exactly; original run metadata also matches after removing
only the new export declaration. Scalar calculations verify all 1680 gene effects,
literal sample/gene axes and reference provenance. Both modes export identical
direct gene effects. Replay 23005042 completes in 9s with 209756 KiB process RSS:
432 tumour profiles across nine frozen refits and 11453 genes reproduce exactly.
These are 48 exposed patients with repeated reference refits, not 432 independent
biological specimens. No new DNA score, comparator selection or claim of improved
validation performance is introduced.

All four jobs exit zero on compute nodes. Source, environment, wheel identity,
tests, parent/child manifests, old/new outputs and replay hashes are retained
under `dev/bulk_gene_export_20261108/{baseline,candidate_v1,probe_source_v1,
validation_23004495,parity_v1,replay_v1}`. The noneditable wheel SHA256 is
`128d9e79cf961fb45153efc7ea21387337664d51411135643dc4c0dfe39b9bec`.
The source-preservation inventory verifies all 512 baseline files and all 41
preexisting test files unchanged; only bulk.py changes among prior production
Python/Rust files. The record is `dev/bulk_gene_export_20261108/preservation.json`.

Retain the information-preserving export. Scientific promotion stops at that
scope: no new genomic inference default, smoothing method, confidence interval
or baseline-superiority claim is supported. Frozen evidence and historical
decisions are unchanged. Original thinning remains 2/5 with three unresolved
failures; the broader clone error/power gate remains failed and certified
biological null n remains zero. Both the profile-accuracy and clone-reliability
goals remain incomplete.
