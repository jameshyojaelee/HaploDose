# Migration to 0.2.0 / schema 2

Schema-2 cellular runs now record final pairwise depth-balance diagnostics
and tree availability under `qc.json`'s `pairwise_depth_diagnostics`.
An untestable comparison withholds that specimen's average-linkage tree.
Previously its missing distance could enter the tree as zero. The optional
`weighted_log_depth_wasserstein` value measures residual log-depth imbalance
after coarse matching. Split and merge criteria retain their declared rules.
Historical trees remain available in their original result directories.
[Review response and validation](SCIENTIFIC_REVIEW_RESPONSE_2026-11-10.md).

The opt-in `stacked_rna_profile` and `stacked_count_profile` APIs now accept
`average_segments=True` together with `gene_segments=True`. This averages
gene-level segmentations under the existing working loss before stacking.
Nonboolean flags and averaging without gene segments fail explicitly. Previous
defaults, profiles and diagnostics reproduce exactly on all181 regression
cases;both original count-to-profile examples also reproduce exactly. These
remain descriptive natural-log relative RNA profiles, with unvalidated
calibration. Evidence covers pooled-query profiles;per-cell DNA,clone and
bulk/spatial benefits remain unestablished. No result schema or historical
artifact is rewritten. [Implementation and scientific limits](NATIVE_PARTITION_AVERAGE_2026-11-10.md).

Historical wheels, source snapshots and result directories remain the authority for reproducing
historical analyses. The new API deliberately rejects undeclared or transformed measurement scales.
Do not relabel reconstructed or normalized values as raw counts to bypass this check.

Default ribosomal exclusion now uses biological symbols. With Ensembl/custom
count IDs, provide an explicit `gene_symbol` column in the positions table for
eligible features. Missing,blank or placeholder symbols fail explicitly; the
previous symbol regex applied directly to non-symbol IDs did not establish
biological exclusion. Symbols annotate categories without changing count IDs
or aggregating counts. Explicit exclude lists still use the declared count-ID
namespace. Single-cell `--keep-ribosomal` (API `exclude_ribosomal=False`) can
disable category exclusion; bulk default exclusion requires annotation.
Historical input/result reproduction uses the preserved historical wheels.

New expression-only, nonspatial fits coalesce exactly equal finite tumour RNA
profiles while retaining the normal component and aggregate mixture mass.
Historical component counts remain reproducible with their frozen wheels.
This uses no similarity threshold and does not equate RNA profiles with genomic
clones; allele and spatial modes retain their existing merge behavior.
[Correction and validation](EXACT_PROFILE_MERGE_2026-11-09.md).

Clone refinement now applies the existing `min_cells` responsibility-mass
minimum to its returned state as well as intermediate iterations. A last
profile update can otherwise leave an ineligible component in the output.
Such components are removed,retained probabilities are renormalized,and all
returned pools are recomputed without refitting allele or nuisance parameters.
Clone0 remains available. This changes affected component counts;reproduce old
outputs with their frozen wheels. It neither merges similar profiles using a
new tolerance nor establishes genomic identity or calibrated membership.
[Counterexamples and validation](FINAL_MEMBERSHIP_SUPPORT_2026-11-10.md).

When clone refinement removes components, new shared-offset residuals weight
observed counts by the same retained responsibility mass as expected exposure.
Historical runs retain the old whole-specimen numerator through their frozen
wheels. This accounting repair retains the legacy clone-only mean; a broader
conditional nuisance-mean replacement failed recovery and is not a default.
It does not establish calibrated clone confidence or genomic validity.
[Evidence and limits](OFFSET_EXPECTATION_2026-11-09.md).

New profile runs also retain `pooled_count`, `expected_count`, `log_ratio_centre`
and `count_units` per bin, with their conditional RNA interpretation in QC.
They preserve the actual pools underlying `log_ratio`; weighted pools need not
be integers. Historical readers/columns remain; historical files are not backfilled.
Aggregate counts and expectations before applying one half-count offset, rather
than treating a mean clone log effect as a pooled-count log ratio. This does not
establish calibrated uncertainty or better DNA accuracy. [Validation and limits](PROFILE_POOL_EXPORT_2026-11-09.md).

New `run` and `run_bulk` calls require `counts_kind='raw'`, `genome_build`, `coordinate_system`
and `gene_namespace`. The CLI uses the corresponding hyphenated flags. Currently supported coordinate
declarations are `one-based-closed` and `zero-based-half-open`; feature namespaces are `symbol`,
`ensembl` and `custom`. Bundled centromeres are hg38; another build requires explicit centromeres.
Duplicate or missing identifiers, invalid count values and incompatible interval declarations fail.
Bulk sample-column headers are checked before parsing can rename duplicates;
ambiguous or empty original headers now fail explicitly instead of being silently
renamed or excluded during sample selection.

An explicit `AlleleParams.cell_kappa` now applies consistently to cell membership,
allele split proposals and balanced-cell candidate selection. The latter two
previously ignored it and used binomial contrasts. `None` retains the legacy
binomial path. Cell-level concentration must be independently justified; pooled
clone concentration is a different quantity. This correction does not establish
biological calibration or support an allele-enabled expression-only null filter.
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
Bulk median mode transforms raw gene counts into derived pseudocounts. Its
segment `count_units` and sample `binned_count_units` now say
`median_derived_pseudocounts`; these totals are not input measurement totals.
The historical `binned_counts` values are retained. Bulk samples additionally
export `input_binned_counts`, the actual sum over retained input genes, with
`input_binned_count_units='input_count_units'`. Sum mode retains input-count units.
Bulk output additionally preserves direct, unsmoothed gene RNA effects in
`gene_rna_log_ratio.npy` (samples by genes, in `samples.tsv`/`genes.tsv.gz` order).
`gene_ratio_metadata.tsv` records mapped input/reference totals and centering.
These float64 effects use pooled normal raw gene rates and a half-count offset,
centred on the median of within-bin gene-effect medians. Each normal excludes
itself from the reference, conditional on the common all-normal feature universe.
No programme effect, HMM dosage state or calibrated uncertainty is transferred.
This retains information lost under intermediate median aggregation; it does
not establish superiority over the naive baseline or DNA copy-number validity.
Bulk runs may request `rna_profile_smoothing='predictive'` (CLI
`--rna-profile-smoothing predictive`). This adds `rna_log_ratio_predictive.npy`
in samples-by-bins order and `rna_profile_smoothing.json` with conditional
masked-gene selection diagnostics. The default is `none`. The descriptive
natural-log RNA effects do not replace HMM calls or acquire calibrated intervals.
Selection conditions on upstream features, normalization and reference estimates;
it is not independent confirmation. Arm boundaries and genomic gaps reset the
smoother. All-singleton inputs retain their bin observations and abstain from
selecting a predictive width. Independent biological benefit remains unvalidated.

Sum-mode bulk runs may instead request `rna_profile_smoothing='native-regions'`
(CLI `--rna-profile-smoothing native-regions`). This allocates fitted bin
expectations to genes using their within-bin reference count fractions, centres
the resulting gene RNA effects, and takes gene medians across contiguous native
RNA-HMM regions. Arm boundaries and gaps reset regions. It adds
`rna_log_ratio_native_regions.npy`, `gene_rna_log_ratio_programme.npy`,
`gene_fitted_expected.npy` and `rna_native_regions.json`, with axes and units in
`run.json`. Median pseudocount mode is rejected because the fitted expectation
does not have the required gene count units. Defaults and existing outputs stay
unchanged. Boundaries and nuisance estimates use the full RNA observations;
these descriptive profiles supply neither independent confirmation nor calibrated
intervals. Development comparisons support a correlation gain over the retained
gene predictor; superiority over strong pseudobulk remains unresolved.

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
Confirmation abstains when multiple candidates share the exact maximum proposal
score. Candidate array order no longer resolves these ties. This is an exact
computational rule, not a calibrated near-tie threshold. Cell exports add
`confirmation_unresolved_reason`, preserving `tied_candidate_scores`,
`no_proposal_depth` or `no_finite_candidate_score`; unprocessed confirmation cells
retain `not_fitted`. Existing `confirmation_status` categories are unchanged.

Candidate groups, allelic imbalance and similarity trees retain exploratory interpretations.
The research whole-site allele splitter enforces tested computational isolation,
but is not integrated into production inference. Reserved SNP evidence does not
by itself validate genomic groups: the [shared-state challenge](ALLELE_SHARED_STATE_2026-11-08.md)
produces repeated reserved contrasts from one DNA population. Independent reads
can retain dependence through allele-expression state. No new confirmation test
or joint-null support is enabled by those experiments.
Absolute copy number, purity, genome doubling and evolutionary ancestry are not identified by
renaming RNA-model parameters. Independent evaluation requires new eligible specimens; former
evaluation cohorts remain regression-only. Read the scientific contract before making a biological
or calibration claim from schema-2 output.

The 0.2.0 candidate repairs the initial thinning-mode counterexample. The later
five-case regression panel remains 2/5, with three unresolved failures; the
joint/full-membership research candidate fails the fresh false-split/power gate
and is not promoted. Interval-coverage failures also remain. Versioning the
output does not make these claims validated. See the
[follow-up evidence](RECOVERY_2026-11-07.md).
Retain unsegmented ratios and per-bin posterior effects alongside segment summaries: averaging
can remove supported focal signal while improving genome-wide correlation. The evidence report
records the exact development cases, operator controls, confidence families and unresolved gaps.

Spatial graph correction (2026-11-09): the public pipeline now canonicalizes
neighbour searches by explicit cell ID. Tied-coordinate reruns may change because
prior ties depended on row order; preserve old outputs and use a new result
identity. `neighbour_graph(...,cell_ids=...)` validates IDs before indexing;
helper calls without IDs retain legacy behavior. This is computational identity
repair, not new clone calibration. All existing spatial assertions and the full
suite pass unchanged. [Evidence](SPATIAL_IDENTITY_2026-11-09.md).

Bin identity correction (2026-11-09): genes with the same genomic start now use
canonical gene ID to determine their order, rather than count-column index.
Reordering equivalent count columns therefore preserves named bin membership and
regional counts. Tied historical bins can change; retain their frozen outputs
and give corrected analyses new result identities. Unique-start controls remain
exactly unchanged. This repairs computational identity without changing bin sizes,
scientific criteria or biological interpretation.
[Counterexample and validation](BIN_IDENTITY_2026-11-09.md).

Bulk per-bin evidence (2026-11-09): new runs additionally retain
`rna_state_posterior.npy`, `rna_log_effect.npy` and `rna_state_index.npy`, with
literal axes and interpretation in `run.json`. These preserve actual conditional
RNA-HMM evidence before segment averaging; they are not biological DNA confidence
or calibrated uncertainty. Existing artifacts remain numerically identical.
Historical runs are not backfilled. [Exact replay and scope](BULK_STATE_EVIDENCE_2026-11-09.md).
