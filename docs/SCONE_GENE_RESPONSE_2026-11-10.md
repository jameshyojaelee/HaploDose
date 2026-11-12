# Gene noise improves the mean profile, but does not pass continuation

[Observed pattern, not a guarantee] A conditional gene NB response raises mean
whole/cell-specific Pearson on the exposed scONE specimen from 0.573712/0.297868
to 0.602558/0.308225. It improves whole profiles in 370/390 cells against the
preceding anchored candidate. The prespecified continuation gate nevertheless
fails: whole-profile wins against pseudobulk are only 181/390, and reference-refined
REFCON still leads both means. Preserve this mechanism result without promoting
a default, expanding its study under the failed gate or changing that gate.

## Scientific question and matched analysis

Regional sums give abundant genes substantial influence. A highly expressed gene
with variable within-region expression can therefore distort the inferred regional
response. This experiment tests whether RNA-estimated gene noise can reduce that
influence beyond an otherwise matched Poisson response. It does not assume that
gene noise distinguishes genomic alteration from coordinated transcription.

The estimand is a cell-specific regional multiplicative RNA response relative
to an other-fold query gene profile, combined with the previously fixed shared
gene-trunk reference contrast. Conditional gene NB scores use means
`N_i * mu_b * q_g * exp(delta_ib)`. Gene dispersions come from variable-depth
fraction moments within the region in the training cells, removing regional total
variation before estimating within-region composition noise. This is a working
moment model: conditionally composed genes are not jointly independent NB.

The log-response score adds `a*(1-exp(delta))`, from the log-response density of
a Gamma(shape=a, rate=a) response prior. This targets a mode in log-response
units; in the Poisson submodel its exponential equals the conjugate posterior
mean response. It is not the response-space MAP. Two prior strengths, one-half
and the training regional moment shape, cross with zero versus fitted gene noise.
All four arms have the same cells, genes, RNA folds, reference, trunk, regional
axes, five-region within-block median and scorer.

Identifying assumptions include transferable gene composition, useful noise
moments and a common regional RNA response. The strongest alternatives are
coordinated expression state, reference mismatch and gene-specific dosage response.
Those can remain CNA-shaped or cause true dosage-responsive genes to be downweighted.
Independent RNA replication and corresponding DNA at the intended resolution
would distinguish those explanations. Optimizer convergence cannot do so.

The fixed specimen has 390 tumour cells, 450 annotated reference cells and 572
5 Mb/arm/gap regions. Smoothing can span five regions and does not establish
5 Mb event detection. RNA is rounded deposited expression, not verified raw
counts. DNA is processed corresponding-cell relative depth, an uncertain
measurement rather than error-free copy states or clone labels. Biological
replication remains one specimen. Earlier evaluation exposure is not erased
by locking this execution's RNA profiles before its DNA scoring.

## Every new arm and the necessary comparators

| Method | Whole Pearson | Cell-specific Pearson | Whole cell wins vs pseudobulk |
|---|---:|---:|---:|
| Existing HaploDose stack | 0.419387 | 0.226796 | See unchanged original table |
| Pseudobulk | 0.561258 | -0.048242 | Reference |
| Previous anchored expected log abundance | 0.573712 | 0.297868 | 161/390 |
| Poisson, half-shape | 0.550747 | 0.302246 | 147/390 |
| Gene NB, half-shape | 0.569637 | 0.317251 | 158/390 |
| Poisson, moment-shape | 0.572777 | 0.297900 | 159/390 |
| Gene NB, moment-shape | 0.602558 | 0.308225 | 181/390 |
| Released REFCON | 0.432251 | 0.371134 | See unchanged original table |
| Reference-refined REFCON | 0.615822 | 0.376007 | See unchanged original table |

The complete report retains 24 methods, all 390 cells, all four numerical score
fields, and 11,260,470 common gene/cell pairs. It does not drop a weaker arm or
change the scoring population. All 7,800 preceding cell records match exactly.
All 96 method/metric means are independently recounted with `math.fsum`.

For the moment-shape gene NB arm:

- Versus its matched Poisson arm: mean whole/cell-specific gains
  +0.029781/+0.010325, with 371/390 and 280/390 cell wins.
- Versus the previous anchored candidate: +0.028846/+0.010357, with 370/390
  and 277/390 wins. Whole-profile median paired gain is +0.027682.
- Versus pseudobulk: +0.041300/+0.356466, with 181/390 and 377/390 wins.
  The whole-profile median paired effect is **-0.006884**. A higher mean
  therefore does not establish improvement for the typical cell.
- Versus reference-refined REFCON: -0.013264/-0.067783, with 166/390 and
  52/390 wins. The competitor gap remains.

The half-shape gene NB arm has the best cell-specific mean among the new arms,
but its whole-profile mean declines relative to the preceding anchored candidate.
Do not combine different arms' best endpoints into a single method claim.
No biological-unit p-values or superiority intervals are supplied for this
single-specimen development comparison.

## Numerical correctness, missing support and scientific decision

The root score is strictly decreasing. Analytically derived global response
brackets, safeguarded Newton steps and binary64 bisection avoid a fitted effect
cutoff. Exact computed zeros or adjacent representable bracket endpoints stop
the solve. Operation/magnitude-based arithmetic diagnostics are saved; these
are not rigorous libm enclosures or biological confidence intervals.

Twelve scalar roots agree with independent exact-rational bisection in response
units. Neutral response and row/gene permutations pass. Six deterministically
RNA-selected real-data roots agree with separate 80- and 160-digit decimal
bisection; maximum discrepancy is 5.55e-17 and the decimal refinements agree
at binary64 precision. The largest saved full-study log-response numerical
diagnostic is 8.48e-12. Independent variable-depth moment fixtures also pass.

Four of 1,716 fold/region noise fits lack two training cells with regional counts
above one; the locked fallback records absent gene-noise evidence and uses the
Poisson submodel. Four moment-prior fold/region fits have zero estimated regional
variance and explicitly pool. These are disclosed lack-of-support boundaries,
not evidence of genomic homogeneity or calibrated certainty.

Continuation required both mean endpoints above the preceding candidate AND
majority whole-profile wins against pseudobulk. All four arms fail. The gene
noise ablation supplies positive evidence for its conditional mean benefit,
but this experiment does not pass its continuation or REFCON-superiority gate.
No larger gene-response study, public default or recovery transfer is supported
by a passing result here. Preserve all negative decisions and do not tune widths,
weights, seeds, exclusions or thresholds to repair them.

The original clone recovery, continuous-state false-split, power and selected-
interval uncertainty problems remain separate. No genomic probability,
specimen error guarantee or interval calibration follows from this profile gain.

## Compute and preservation

Both jobs use nslab account/QOS, CPU partition, one CPU, existing Conda Python,
new output directories and v4 Intel placement matching previous scoring.
Package/native files, production dependencies and original tests are unchanged.
Full package/Rust/GPU suites were not run for these research-only additions.

| Job | Work | Progress | RAM / time requested | State and elapsed | Measured process peak RSS |
|---|---|---|---|---|---:|
| 23192397 | Numerical checks, four arms, fixed DNA scoring | 3/3 folds ████████ | 8 GB / 1 h | Completed 0:0, 2m48s | 928,352 KiB |
| 23193489 | Independent roots and paired-cell decisions | 6/6 roots ████████ | 4 GB / 15 min | Completed 0:0, 26 s | 388,672 KiB |

Two completed, zero failed, zero running or pending in this experiment. The
separate clone preflight 23191239 is left untouched; consult its own journal
for live state. Scheduler memory sampling differs from the process peaks and
is retained separately. Earlier profile and failed clone artifacts are not changed.

[Prespecified design](../dev/scone_gene_response_20261110/DESIGN.md),
[root solver and rational references](../dev/scone_gene_response_20261110/response.py),
[executed fitting/scoring](../dev/scone_gene_response_20261110/study.py),
[independent reader](../dev/scone_gene_response_20261110/summarize.py),
[source hashes](../dev/scone_gene_response_20261110/source.sha256),
[reader hashes](../dev/scone_gene_response_20261110/reader_source.sha256),
[RNA predictions locked before DNA](../dev/scone_gene_response_20261110/run_23192397/RNA_PREDICTIONS_LOCKED.json),
[24-method score report](../dev/scone_gene_response_20261110/run_23192397/REPORT.json),
[all cell scores](../dev/scone_gene_response_20261110/run_23192397/scored_cells.tsv),
[preceding records exact](../dev/scone_gene_response_20261110/run_23192397/PREVIOUS_EXACT.json),
[independent checks and all decisions](../dev/scone_gene_response_20261110/summary_23193489/results/SUMMARY.json),
[all paired comparisons](../dev/scone_gene_response_20261110/summary_23193489/results/paired_cells.tsv),
[fit submission](../dev/scone_gene_response_20261110/submission.json),
[reader submission](../dev/scone_gene_response_20261110/reader_submission.json),
[full scheduler accounting](../dev/scone_gene_response_20261110/accounting_complete.tsv),
[prior ledger entries preserved](../dev/scone_gene_response_20261110/ledger_preservation.json),
[preceding candidate](SCONE_REGIONAL_SHRINK_2026-11-10.md).
