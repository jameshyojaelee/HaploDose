# Default ribosomal exclusion follows biological symbols across count-ID namespaces

COMPLETE:all3jobs exit0:0. The candidate passes493Python and4Rust tests,
with1GPU-only skip and8existing warnings. Both old-wheel equivalence failures
are corrected; all10new annotation-contract tests pass. Existing assertions
and numerical criteria remain unchanged. Environment preparation takes4m13s;
candidate build/tests3m27s,including164.22s pytest. Native binary is unchanged.
[Python tests](../dev/ribosomal_namespace_20261110/validation_23160780/pytest.txt),
[Rust tests](../dev/ribosomal_namespace_20261110/validation_23160780/cargo_test.txt),
[installed artifact](../dev/ribosomal_namespace_20261110/validation_23160780/installed_artifact.json),
[accounting](../dev/ribosomal_namespace_20261110/accounting_completed.tsv),
[preservation](../dev/ribosomal_namespace_20261110/preservation_completed.json).

Earlier preparation checkpoints retain their original context.

Baseline23160778 completes0:0 in13s:the inner pytest fails both intended
equivalence assertions,with zero execution errors. This confirms the existing
namespace defect on the installed wheel. Candidate build/full-suite validation
awaits the live environment job; no passing candidate result yet.
[Old-wheel counterexamples](../dev/ribosomal_namespace_20261110/baseline_23160778/report.json),
[full assertion details](../dev/ribosomal_namespace_20261110/baseline_23160778/pytest.xml).

Implemented review finding13 as an explicit annotation contract. Previously
the symbol regex was applied directly to all count IDs. Ensembl IDs could not
match ribosomal symbols,so biologically identical inputs selected different
features. No inference run is needed to establish that source behavior.

Now default ribosomal exclusion uses count IDs directly only in the symbol
namespace. Ensembl/custom IDs require an explicit `gene_symbol` positions
column for eligible features. Missing,blank,placeholder or nonstring required
annotations fail explicitly. Annotation does not change count IDs,strip ID
versions,guess aliases or merge counts. A custom ID resembling an RP symbol
is retained when its declared biological annotation is nonribosomal. The
used-gene output retains the symbol mapping for this category-filter path.
Other explicitly supplied regexes and exclusion lists keep operating on count
IDs. Explicitly disabling category filtering remains possible.

This is a versioned input migration,documented in README,CLI help and
MIGRATION_0.2.md. It may change non-symbol feature selection or reject inputs
that previously ran with ineffective biological filtering. Frozen code and
results remain authoritative for their historical analyses. Symbol-namespace
filtering is unchanged. No new scientific threshold or biological performance
claim is introduced.

The new tests require exact retained-gene and bin-mapping equivalence for
symbol,Ensembl and custom representations. They also check invalid/missing
annotations,explicit bypass and unused annotations. Existing synthetic
custom-ID fixtures receive only their declared synthetic symbol column;
count values and existing assertions are unchanged. AST comparison confirms
every existing assertion is identical. No numerical tolerance or recovery
criterion changes.

Submitted three CPU jobs:no arrays or interactive QOS.

| Job | Purpose | Account/QOS | CPU | RAM | Wall |
|---|---|---|---:|---:|---|
|23160778|old-wheel equivalence counterexamples|nslab|1|4GB|30min|
|23160779|private conda environment|innovation|1|8GB|1h|
|23160780|candidate wheel and full Python/Rust tests|nslab|2|32GB|1h|

Validation depends afterok on environment preparation. Baseline runs against
the previously validated allele-corrected wheel; the new counterexample tests
are read from the immutable source_before. Both snapshots contain924files.
The original error/power worker23159097 and dependent summary23159099 continue
unchanged. Executable results are pending; syntax and assertion preservation
are static checks,not passing inference or scientific-validation results.

[Binning correction](../python/haplodose/bins.py),
[new tests](../tests/test_ribosomal_namespace.py),
[migration](MIGRATION_0.2.md),
[locked design](../dev/ribosomal_namespace_20261110/DESIGN.md),
[source differences and assertion preservation](../dev/ribosomal_namespace_20261110/source_difference.json),
[all commands and handles](../dev/ribosomal_namespace_20261110/jobs.json).
