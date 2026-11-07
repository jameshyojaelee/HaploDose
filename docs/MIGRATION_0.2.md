# Migration to 0.2.0 / schema 2

Historical wheels, source snapshots and result directories remain the authority for reproducing
historical analyses. The new API deliberately rejects undeclared or transformed measurement scales.
Do not relabel reconstructed or normalized values as raw counts to bypass this check.

New `run` and `run_bulk` calls require `counts_kind='raw'`, `genome_build`, `coordinate_system`
and `gene_namespace`. The CLI uses the corresponding hyphenated flags. Currently supported coordinate
declarations are `one-based-closed` and `zero-based-half-open`; feature namespaces are `symbol`,
`ensembl` and `custom`. Bundled centromeres are hg38; another build requires explicit centromeres.
Duplicate or missing identifiers, invalid count values and incompatible interval declarations fail.
These declarations record provenance; software cannot establish that a supplied declaration is true.
Allele-enabled runs additionally require `allele_genome_build` (`--allele-genome-build`) to match
the expression annotation assembly. A copied VCF header alone is not assembly verification; trace
the aligned reference and executed variant-panel, map and annotation inputs before declaring it.

Schema-2 runs write `run_manifest.json` linking input hashes, parameters, installed source/native
binary, environment and outputs. Failed execution is recorded when failure occurs after run creation.
Multi-dataset inputs yield independent subruns indexed by `datasets.json`; they no longer share a
feature universe or nuisance estimates. A parent result directory therefore need not contain a
single combined cell/profile table. Use each child's manifest and QC.

Internal and exported intervals are zero-based half-open. Numbat SNP `POS` remains one-based and is
converted at interval comparisons. Segment tables retain historical numeric columns for readers,
but schema-2 effect and uncertainty columns state RNA log-ratio units explicitly. Historical names
do not grant their former biological interpretation. Approximate standard errors exclude adaptive
selection and fitted nuisance uncertainty. Arm-score diagnostics are unvalidated; the misleading
`arm_calls_calibrated` flag and `min_detectable_log_dose` field are removed. The corresponding
new scale diagnostic is `scaled_se_rna_log_ratio_approx`, with `uncertainty_scale_z` recorded;
it is not a detection limit. Existing readers can still read unchanged historical result files.
Use `pooled_count` and `expected_count` for measurement totals. Their units are the source count
units, which may be UMIs, reads or fragments; the deprecated `pooled_UMI`/`expected_UMI` aliases
do not establish that reads have been deduplicated to molecules.
Cell tables similarly provide `library_count` and `binned_count` alongside the historical names.
`observed_mean_rna_log_ratio` is the centre corresponding to the approximate log-ratio SE.
It is distinct from `rna_log_effect`, the HMM posterior mean state effect. The export records
`uncertainty_statistic` and `uncertainty_center_field`; `uncertainty_estimand` is null
because the approximate scale diagnostic is not a validated parameter interval.
This corrects the earlier naming of an observed statistic as an estimand and does
not validate an interval around the posterior effect. Historical files are unchanged.

Reference-only feature selection and independent dataset processing can change outputs even with
the same seed. Query/reference mismatch is now recorded as a diagnostic rather than used to exclude
possible signal. Allele-balanced query cells are no longer automatically admitted as references;
the explicit `--in-sample-reference` option is exploratory. Allele-enabled inference combined with
the expression-only split null is rejected because the observed and simulated statistics differ.
Gene-veto-enabled runs also reject this filter because its bin-only simulations omit the veto.
The public gene-aware pipeline can still run with the default `split_null_draws=0`; it does not
claim specimen error control. Low-level bin-only filter experiments are diagnostic-only.

Optional query roles can be declared by metadata `analysis_split` values `discovery` and
`confirmation`, or allocated with `--discovery-fraction` and the fixed seed. Use one mechanism.
References retain their own fit/holdout roles. Confirmation values cannot modify discovery fits,
features or proposals. Frozen expression-only candidate scoring and signed contrast diagnostics
are exported separately. Disjoint feature sets do not prove independence of correlated biological
measurements, so these diagnostics are not calibrated confirmation tests or validated DNA membership.

Candidate groups, allelic imbalance and similarity trees retain exploratory interpretations.
Absolute copy number, purity, genome doubling and evolutionary ancestry are not identified by
renaming RNA-model parameters. Independent evaluation requires new eligible specimens; former
evaluation cohorts remain regression-only. Read the scientific contract before making a biological
or calibration claim from schema-2 output.

The 0.2.0 candidate repairs the original thinning-mode recovery regression but retains
failures on other fixed simulation draws and broader clone-error/interval-coverage
experiments. Versioning the output does not make these claims validated. See the
[follow-up evidence](RECOVERY_2026-11-07.md).
Retain unsegmented ratios and per-bin posterior effects alongside segment summaries: averaging
can remove supported focal signal while improving genome-wide correlation. The evidence report
records the exact development cases, operator controls, confidence families and unresolved gaps.
