# Kinker measurement chain and complete-population scoring

The two locked development cases conserve counts and have complete exported
profile support. Missing-profile reweighting therefore does not explain their
profile gap. A separate reproducible scorer defect is corrected: unavailable
clones or bin values can no longer become a smaller scored population or zero.
No inference default, frozen result or scientific criterion changed.

## Primary measurement chain

The cases were already fixed in `benchmarks/hardening/claims.json`:
`2313287_STOMACH` and `42MGBA_CENTRAL_NERVOUS_SYSTEM`. Their biological units are
two cell lines, not their 620 query cells or the number of checked matrix entries.
They are mechanism probes, not an independent performance panel.

| Case | Rotation | Query cells | Fit-reference cells | Genes / bins | Exact count entries |
|---|---:|---:|---:|---:|---:|
| 2313287_STOMACH | 0 | 103 | 2118 | 10133 / 1003 | 2227663 |
| 42MGBA_CENTRAL_NERVOUS_SYSTEM | 1 | 517 | 2467 | 10116 / 1001 | 2986984 |

Input is `results/kinker_input_v2_20261101_1551/counts.h5`; outputs are
`results/kinker_field_20261102_0529/haplodose_r{0,1}`. New evidence lives under
`dev/kinker_measurement_chain_20261109/`. `chain_v1/identities.json` records exact
ordered query/reference cell IDs, matrix rows, genes, columns and gene-to-bin
mapping. Its manifest hashes counts, truth, positions, metadata, exported profiles,
bin counts, QC and the inspected preparation/scorer source. `support_v3` verifies
all input and output hashes from that audit again.

Independent sparse aggregation exactly reproduces all 5,214,647 selected bin
entries. Full matrix row sums exactly reproduce metadata UMI values. All 103
and 517 query cells are assigned; each required profile value is finite. The
first line has one exported candidate group, the second six. Cell-line origin
does not establish one genomic clone or validate these six groups.

The largest absolute difference between full-library and retained-gene reference
line weights is 0.0058323552 and 0.0101816208, respectively. These are differences
in compositional weights, not measured errors. A DNA-favouring reweighting is
not justified. The current historical DNA contrast uses metadata UMI-weighted
DepMap values; its relationship to the pooled RNA reference still requires a
declared observation model and verified DNA measurement scale.

## Correction and independent falsification

`benchmarks/kinker/evaluate.py` is revision `1.2_complete_population`.
Previously, missing clone profiles were excluded and remaining clone weights
renormalized. Pandas summation could also turn an all-missing value into zero.
The corrected scorer requires every query clone and scored bin to have a finite
value; otherwise that specimen/method is unavailable, with represented fraction
and failure reason exported. An unavailable sign score remains unavailable.
This is a completeness contract, not a new biological threshold.

The estimand remains the cell-weighted arithmetic mean of exported clone log
effects. It is explicitly not a pooled-count log ratio. Comparator tests share
one complete population across the available method family. Undefined tests stay
in the Holm family for adjustment, with no reported p-value or success verdict.
Coverage and complete-family population size remain visible in the summary.

Falsification fixtures remove a clone, a bin, an assignment, or a finite value;
none may yield a numerical profile for the full population. Other fixtures check
independently calculable weighted means, row/bin permutations, duplicate rejection,
shared comparator populations and an entirely unavailable comparator. All original
input-contract tests remain intact. Final compute run: **46 passed**.

For both locked cases, the historical and corrected profile aggregates agree
exactly across three quantities and 6,012 bin-level calculations. Independent
scalar `math.fsum` references differ by at most 2.220446049250313e-16, within a
binary64 operation-count bound based on absolute weighted sums. No arbitrary
universal epsilon or relaxed scientific boundary is used.

## Identification and validation remain separate

- **Computation:** selected count and identity chains pass; the missing-support
  scorer defect has a reproducer, fix and independent regression checks.
- **Identifiability:** counts were reconstructed from published CPM using each
  cell's smallest nonzero value. Integer recovery does not prove original molecule
  counts: missing observed singletons can leave a cell's scale ambiguous. Original
  counts/libraries are needed to settle this. Cell-line assignment was inferred
  from expression similarity to DepMap; source tags or genotypes would independently
  establish correspondence. Retained-feature composition, RNA regulation and DNA
  dosage remain competing explanations for profile differences.
- **Biology:** DepMap is a separate assay/culture. There is no same-cell DNA truth
  here. This audit does not validate CNA membership, confidence coverage or a
  specimen false-split guarantee. Allele/pileup evidence was not investigated:
  no count/barcode discrepancy in this bounded audit triggered read inspection.

The exact historical executed source/environment is not newly established by
hashing today's inspected code or old QC. Current isolated wheel execution is
recorded for the audit; the missing historical execution provenance remains a gap.
The truth scale, specimen correspondence and original molecules require evidence,
not additional simulations. No model correction is proposed from this negative
diagnostic, and no DNA-driven reference selection is performed.

## Compute evidence

All executable checks ran on compute nodes with immutable source snapshots and
the validated noneditable wheel from `dev/bulk_posterior_export_20261109`.
No package/native source changed in this work; the scorer, its tests and audit
scripts changed. Jobs: 23069935 (innovation, 2 CPUs, 16 GB, 13 s), 23070330
(nslab, 2 CPUs, 8 GB, 8 s), 23070516 (nslab, 2 CPUs, 8 GB, 6 s), all exit 0.
The intermediate successful test/audit outputs remain frozen alongside the final
ones. Resource and execution records are in `analysis_JOB/`; `jobs.json` records
requested resources and process peak RSS. No downstream job is pending.

The original three thinning failures, failed false-split/power tradeoff and
uncertainty failures remain unresolved. Certified biological clone-null n remains
zero. These results authorize neither a superiority claim nor weaker criteria.
