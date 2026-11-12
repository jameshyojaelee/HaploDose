# Deliver the improved conditional RNA stack as an opt-in package API

The installed wheel passes431 Python tests and4 Rust tests, with one GPU-only
skip, and reproduces all181 frozen RNA profiles and diagnostics exactly. Field
scoring is complete. No independent biological
validation, clone5% guarantee, bulk improvement or universal superiority is
claimed. [Prior frozen scientific evidence](RNA_STACKED_SMOOTHING_2026-11-10.md)
establishes the task-specific conditional regression gain and failed bulk transfer.

## Implementation and scientific contract

New [package module](../python/haplodose/stacking.py),
[numerical/invariant tests](../tests/test_stacked_smoothing.py),
[installed parity runner](../benchmarks/hardening/stacked_package_parity.py),
[delivery design](../dev/native_stacked_profile_20261110/DESIGN.md).

`stacked_rna_profile` accepts finite prepared natural-log RNA gene contrasts,
complete zero-based integer bin IDs and contiguous arm/gap labels. It validates
the original bin IDs before conversion, so floats and booleans cannot be
silently truncated. Existing median-selector validation and fallback remain.
The proven research expert bank, masking, operation order, sparse L1 objective,
feasible-dual check and numerically unresolved fallback are retained.

Tests independently enumerate exact rational L1 breakpoints and direct-residual
segment partitions, verify target masking and genomic resets, check constant
fallback/singleton abstention, and reject malformed identities. Existing tests
and all scientific criteria are unchanged. NumPy/Pandas/SciPy are existing
production dependencies; no dependency was added.

Estimator: regional descriptive RNA abundance relative to an explicit reference.
Working Gaussian segments are smoothing experts, not raw-count likelihoods.
Target values train the mixture weights; masked base predictions do not create
independent confirmation. Upstream fitting is conditional. Diagnostics identify
natural-log RNA units and unvalidated calibration. RNA-equivalent regional
transcription/reference effects remain competing explanations; independently
corresponding DNA is required for a genomic claim.

This is an opt-in API rather than an automatic change to bulk or clone decisions.
The evidence for this change is improved single-cell profile prediction, while
the same fixed stack worsens bulk performance. Current confirmed/merged clone
rules, original recovery assertions and5% specimen confidence target are intact.

## Compute chain and reproducibility

Seven CPU/v6&amd jobs were submitted, no arrays. All writes are to the new task,
and all analysis starts from a neutral directory with recorded artifacts.

| Job | Role | Account/QOS | CPUs/RAM/time | Dependency |
|---|---|---|---|---|
|23154878|Exact isolated dependencies, no shared install|innovation|1/8GB/1h|none|
|23154879|Noneditable wheel, full Rust/Python tests|nslab|4/32GB/2h|afterok23154878; pass|
|23154880|Original replay launch|innovation|1/8GB/2h|exit2 before analysis|
|23154887|Original field launch|innovation|1/16GB/2h|exit2 before analysis|
|23154909|Corrected launcher,181-line installed replay|innovation|1/8GB/2h|complete, exact|
|23154910|Field identity-contract failure|innovation|1/16GB/2h|exit1 before RNA/DNA scores|
|23155009|Corrected population-contract field comparison|innovation|1/16GB/2h|complete|

[Exact commands](../dev/native_stacked_profile_20261110/jobs.json),
[complete source_v1](../dev/native_stacked_profile_20261110/source_v1/source_manifest.json),
[field design frozen before scores](../dev/native_stacked_profile_20261110/FIELD_DESIGN.md).
Source_v2 adds the field runner. Source_v3 corrects its population-record
interpretation; package code/tests are unchanged fromv1. The test job originally
requested8 CPUs, then all scheduler CPU/task/node fields were reduced to4 on
the identical AMD hardware to fit available capacity; RAM/time/tests unchanged.

Both original follow-ups exited2 because the previous launcher's exact allowed
environment path still named the preceding validation environment. A new task-
local copy changes only that path; wheel/source/import/protected-path guards
remain. [Preserved correction](../dev/native_stacked_profile_20261110/LAUNCHER_CORRECTION.md).
The next field execution finds an actual scorer-contract mistake: the frozen
`ordered_cells` list is query PLUS fit-reference. The corrected scorer reconstructs
that full ordered subset, then requires exact ordered query/ref populations and
reference DNA weights across every method. No field outcome was available before
the fix. [Exact source evidence and failed record](../dev/native_stacked_profile_20261110/FIELD_IDENTITY_CORRECTION.md).

The field comparison keeps the five original RNA baselines, gene-level
pseudobulk, frozen HaploDose gene effects, inferCNV R, infercnvpy, CopyKAT and
SCEVAN. Every line uses identical query cells and common gene/region support.
Complete-population profile averaging prohibits missing groups from improving
a score by disappearing. Two endpoints,1MbPearson and armSpearman, form a fixed
family22 with20k shared-reference-group bootstrap draws, seed2026111002.
Simultaneous basic intervals and Holm p<.05 must both pass. An incomplete
comparator family prevents a superiority verdict. Fixed upstream fitting and
eight dependence groups provide only approximate conditional uncertainty on
exposed regression; this is not new independent evaluation or calibrated coverage.

## Outcomes

Full [Python suite](../dev/native_stacked_profile_20261110/validation_23154879/pytest.txt):
431 passed,1 GPU-only skipped, two visible HiGHS option pass-through warnings,
147 seconds. [Rust suite](../dev/native_stacked_profile_20261110/validation_23154879/cargo_test.txt):
4 passed,0 failed. No test/tolerance changed. Test-process peakRSS2,263,236KiB.

[Installed artifact](../dev/native_stacked_profile_20261110/validation_23154879/installed_artifact.json):
noneditable wheel SHA256 `b6d0210c7d3e461c9e320123dcd6027d1189a4bbdbb45cde10cba0be2eb6bf87`;
native SHA256 `e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`,
identical to the preceding native binary. The new Python module is the intended
package addition, not a Rust or clone-model change.

The installed [181-line replay](../dev/native_stacked_profile_20261110/parity_v1/report.json)
passes exact equality for every profile, selected/fallback decision, expert
weight and numerical diagnostic, with no DNA read. Runtime7m33s. No tolerance
was introduced to bridge the package/research result. Every existing79 package,
native and test file is byte-identical to the initial baseline; only the new
module and its new tests extend those directories. All153 earlier ledger fields
remain exactly equal. Public README instructions and evidence docs are updated.

All15 jobs across research and delivery are terminal:12 complete,3 failed
attempts retained and corrected (two launcher paths, one scorer record
interpretation). No failed test was relaxed, no candidate changed after new
field outcomes, and no frozen result was written. Resource/accounting records
and complete commands remain in each task directory. The code is computationally
validated for the tested inputs; default promotion, independent performance,
biological calibration and the broader scientific goal remain unachieved.

## Matched field performance, including unavailable evidence

[Complete181-line score table](../dev/native_stacked_profile_20261110/field_v2/scores.tsv),
[fixed-family report](../dev/native_stacked_profile_20261110/field_v2/report.json),
[population coverage](../dev/native_stacked_profile_20261110/field_v2/population_coverage.tsv).
All query/reference identity, order and UMI DNA-weight checks pass. CopyKAT
has13 incomplete populations/profiles; SCEVAN has14. Their unavailable values
remain in the181-line table and fixed family. Neither the query population nor
the comparator family is reduced to improve a success verdict.

Subsequent bounded inspection checks the first listed incomplete population,
8305C_THYROID, rotation0/pool19. Both native matrix cell lists omit exactly
two of its447 query cells: GCAGTTATCTTGTATC-15-19 and
TGAGAGGTCACTTATC-15-19. CopyKAT's original prediction file explicitly labels
both `not.defined`; both exported calls tables retain missing calls. Thus this
case is not merely a scorer losing otherwise available native matrix columns.
The precise upstream exclusion mechanism has not been traced, and this one
case does not establish the cause of every missing population. No new run,
zero-profile substitution, cell removal or formal superiority verdict follows.

[CopyKAT original predictions](../results/kinker_copykat_20261101_2246/r0_pool19/copykat_20261101_225710/native/run_copykat_prediction.txt)
and [SCEVAN matrix cell identities](../results/kinker_scevan_20261101_2246/r0_pool19/scevan_20261102_002703/native/cells.txt)
preserve the relevant native evidence.

The following are **descriptive complete-family means on167/181 lines**, using
the same all-method genes, windows and arms. They must not be confused with the
all181-line original-support mean0.707226 or older unmatched arm medians.

| Frozen method | Mean1MbPearson | Mean armSpearman |
|---|---:|---:|
| Raw bin pseudobulk |0.594591|0.838167|
| Width-five pseudobulk |0.662165|0.835737|
| MDL pseudobulk |0.664927|0.828644|
| Retained predictive median |0.674665|0.838856|
| Preceding predictive partition |0.681677|0.836294|
| Gene-level pseudobulk |0.351555|**0.849108**|
| Frozen HaploDose gene ratio |0.578710|0.838156|
| inferCNV R |0.668477|0.813294|
| infercnvpy |0.625294|0.798036|
| CopyKAT |0.609403|0.736958|
| SCEVAN |0.568289|0.684473|
| Conditional RNA stack |**0.702970**|0.847446|

At1Mb the stack has the highest point estimate, including139/167 wins against
inferCNV R. At arms, raw gene pseudobulk is slightly better; the stack's paired
mean difference is-0.001661858,82 wins/83 losses/2 ties, simultaneous interval
[-0.008556793,0.005164675]. Do not replace a competitive simpler arm summary
with the fine-scale smoother on the strength of the1Mb result.

The prespecified missing-family rule blocks **every formal field superiority
verdict**:22 tests retain their place, are marked untestable and get adjusted
p=1. No dropping the14 difficult lines/methods, changing families or narrower
post-outcome inference. Positive conditional point differences do not override
this rule. The all181-line five-baseline regression evidence remains a separate
prespecified result; neither comparison provides prospective biological validation.

Matched smoothing can improve descriptive RNA fidelity without creating genomic
information. The RNA-equivalent transcriptional alternative, cross-culture DNA
uncertainty and fixed upstream reference fitting remain. No individual clone
membership or CNA truth is established by these specimen-level correlations.

## Files changed in this performance increment

- Package/API and tests: `python/haplodose/stacking.py`, `tests/test_stacked_smoothing.py`.
- Research implementation and runners: `benchmarks/hardening/rna_stacked_smoothing.py`,
  `benchmarks/hardening/rna_stacked_summary.py`, `benchmarks/hardening/bulk_stacked_smoothing.py`,
  `benchmarks/hardening/stacked_package_parity.py`, `benchmarks/hardening/stacked_field_comparison.py`.
- Public instructions and claim/status records: `README.md`, `benchmarks/hardening/claims.json`,
  `docs/STATUS.md`, `docs/DEVELOPMENT_LOG.md`.
- Detailed reports: `docs/RNA_STACKED_SMOOTHING_2026-11-10.md`,
  `docs/NATIVE_STACKED_PROFILE_2026-11-10.md`.
- Rebuilt navigation: `docs/DEVELOPMENT_INDEX.md`, `docs/DEVELOPMENT_ATTEMPTS.md`,
  `docs/DEVELOPMENT_RUNS.md`. This inventory build does not rerun scientific validation.

New immutable sources, designs, exact job commands, isolated environment/wheel,
all successful/failed logs, predictions, scores, accounting and memory records
are under `dev/rna_stacked_smoothing_20261110` and
`dev/native_stacked_profile_20261110`. Their linked source manifests and
preservation records distinguish this work from all earlier dirty edits.
