# Bulk sample identity and derived-count units

Two remaining input/output correctness defects are corrected. Neither repair
changes the bulk estimation method, its thresholds or the numerical output for
well-formed inputs. The existing scientific recovery and false-call assertions
are retained; independent old-wheel replay verifies exact numerical parity.

## Original sample headers

The previous bulk reader passed the TSV directly to pandas, which can rename
duplicate sample headers before HaploDose selects the requested samples. For
example, two columns both named N0 are ambiguous observations, not evidence that
one should be renamed N0.1 and silently excluded. The corrected reader inspects
the original header first and rejects empty or duplicate sample IDs. It then
requires unchanged header identities after full parsing.

Five counterexamples/invariants cover duplicate, empty and whitespace-only IDs,
and exact preservation of distinct N0, N0.1 and 001 identities in plain/gzip TSV.
No ambiguous column is automatically aggregated or assigned a new identity.

## Observed counts versus derived pseudocounts

Median bulk mode transforms gene log ratios into a derived bin measurement.
Although its run-level likelihood already identified this approximation,
segment totals incorrectly retained `count_units='input_count_units'`, and
sample totals lacked a unit field. These totals need not equal the retained
input observations. In the fixed test's first normal sample, the median-mode
total is7868526.155594692 while the retained raw total is7921955.

Median-mode segments and samples now explicitly use `median_derived_pseudocounts`.
The historical numeric fields remain available. The new sample field
`input_binned_counts` reports the actual sum over retained input genes, with
`input_binned_count_units='input_count_units'`. Sum mode keeps input-count units.
The source assay still determines whether input counts mean reads, fragments or
UMIs; these fields do not infer molecular provenance.

An independent integer sum from the original count table and exported feature
mapping checks the raw totals. The falsifier for a units-only numerical change
is any altered original profile, call or numeric table value for a valid input.
Both sum and median replays pass exact comparisons, without a tolerance.
This establishes computational behavior, not count-model calibration or DNA
validity. The approximate uncertainty and candidate-only interpretations remain.

## Execution and preservation

The noneditable wheel in `dev/bulk_measurement_units_20261108/environment` is
`b93997e10ec40198a35dbccb395aa17bad43936d19b664f44e9c249810103522`.
Environment copy22996512 uses the preceding immutable environment's exact
dependencies and completes in3m37s. Build/validation22996557 completes in7m12s:
330 Python tests pass, one GPU-only test is skipped, and all four Rust tests pass.
Python-test process peak RSS is2505456KiB. All execution is on compute nodes.

Parity job22996562 completes in53s with scheduler peak RSS202000KiB. It replays
both modes with the preceding terminal wheel, verifies five original tables per
mode and both profile arrays exactly, independently checks raw totals, and
verifies446 preserved baseline files. Corrected unit labels and the declared new
metadata fields are the only permitted differences. Package/test changes are
limited to `python/haplodose/bulk.py` and `tests/test_bulk.py`; the source inventory
is recorded in `validated_core_changes.json`.

Candidate v2 is the executed source. Candidate v1 and its audit snapshot remain
preserved; their two pending jobs were cancelled before execution when the header
defect was found. No shared installation or frozen result was overwritten. The
ongoing joint-clone experiment retains its original terminal wheel and source.
See `dev/bulk_measurement_units_20261108/parity_v1/report.json` and the validation
directory for commands, source/environment identities and all check results.
