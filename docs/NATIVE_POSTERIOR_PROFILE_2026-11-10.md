# Retained bulk RNA profile estimator delivered as an optional package API

The scientific review corrected the original within-block runtime claim.
See the [measured scaling and exact replay](SCIENTIFIC_REVIEW_RESPONSE_2026-11-10.md)
for the updated implementation and complexity.

The package now exposes the retained reference-weighted boundary-averaged RNA
profile and an identity-aligned equal-reference mean. Installed outputs reproduce
ALL48 frozen central RNA profiles and fitted decisions EXACTLY; the mean helper
reproduces ALL39 fixed eight-reference ensembles EXACTLY. The isolated wheel
passes477 Python/4 Rust tests,1 GPU-only skip. Existing package modules, original
tests, defaults, clone method and frozen results remain unchanged. This delivers
the preceding performance improvement; it does not create a new superiority,
calibration or biological result.

## What is now usable

Public [posterior.py](../python/haplodose/posterior.py) provides:

- `posterior_rna_profile`:prepared query gene RNA effects/expected exposures,
  corresponding reference effects/exposures, complete gene/bin and contiguous
  arm/gap blocks, explicit reference identity and source measurement declaration.
- `average_rna_profiles`:equal mean of separately prepared profiles, with
  explicit genomic-region identities, fitted-reference identities, alignment and
  identical finite support. Missing fits/regions cannot silently disappear.

The first uses the retained normal-centred reference gene weights, normal
leave-one-out working bin scale, existing bulk beta grid/tau floor, zero-centre
Gaussian boundary model, exact-dyadic interval medians and posterior run weights.
The tested mathematical implementation is ported directly into the package;
research scripts are not imported at runtime. It averages relative RNA medians
over constant-state runs. Its grid and conditional model probabilities are
smoothing devices, not validated DNA states or membership probabilities.

Estimand:regional relative RNA abundance in natural-log units against the
declared normal panel. Observation model:Gaussian working boundaries and
robust gene locations on already prepared RNA effects. Expected gene values
are fitted exposures, not new count observations. Reference effects/expectations
must share gene identity/order and upstream preprocessing with the query.
The caller must establish this correspondence; this prepared-array API cannot
infer missing gene labels. At least3reference profiles are needed for the
existing leave-one-out scales. Reference spread is a working weight, not an
identified technical variance or an interval SE. Quadratic within-block storage
and computation remain; genome-wide dense use without genomic blocks is not
established as scalable.

Identifying assumption:the reference and regional location adequately represent
the relevant RNA contrast. Shared lineage mismatch or coordinated transcription
can explain DNA-shaped RNA effects in every reference fit. Independently
processed corresponding DNA at matching resolution tests genomic prediction;
independent RNA references test shared-reference sensitivity. Cells, regions,
and reference refits are not independent biological replication. This increment
adds0new biological specimens and does not validate DNA states, purity, WGD,
genomic clone membership or spatial assignments.

The ensemble API preserves a crucial nonlinear order:map EACH reference fit
to common intended regions, THEN average. Averaging full-bin profiles before
taking window medians is a different estimator and is not claimed to reproduce
the historical window ensemble. It aligns each region axis to the first, rejects
duplicate/missing identities and support disagreement, and uses compensated
`fsum`. All-NaN common support remains NaN with zero finite regions recorded.
Fit count is not a posterior sample size or biological n.

[Frozen delivery contract](../dev/native_posterior_profile_20261110/PORT_DESIGN.md),
[exact source-kernel identities](../dev/native_posterior_profile_20261110/EXTRACTED_KERNELS.json),
[README usage](../README.md).

## Decisive checks and their outcomes

The delivery criterion was EXACT equality, not a new epsilon. Independent
tests enumerate all27short state paths directly and compare complete-run
probabilities against forward/backward evaluation. An operation-count rounding
bound supplies the numerical allowance. Direct Fraction weighted medians,
midpoint ties and compensated ensemble sums independently check the aggregation.
Genomic-reset factorization, gene/reference/region permutations, reconstructed
label preservation, and invalid/missing/duplicate/support cases pass. No test
mirrors a desired scientific ranking or modifies an old acceptance rule.
[New tests](../tests/test_posterior_profile.py).

The installed public function replays ALL48 central query profiles from the
retained20-normal/68-total RNA contexts. Same ordered genes/bins, expected
exposures, normal membership, genomic resets and ProfileParams. Every profile
and each sample's beta,tau,centre,full grid scores,near-tie decisions and rounding
budget matches exactly. Normal output rows remain unavailable exactly as before.
All48outputs lock before any DNA access; no DNA was read by this delivery job.
[Central report](../dev/native_posterior_profile_20261110/central_v1/report.json),
[locked outputs](../dev/native_posterior_profile_20261110/central_v1/RNA_LOCKED.json).

The independent ensemble delivery job calls the public identity-aligned mean
on all39methods x48specimens x8fixed reference fits at the original window
resolution. Every complete array equals the historical compensated mean exactly,
including unavailable values. The source archives also contain DNA members,
which this job does not open. No selected method, region, patient or reference
draw is dropped to pass.
[Ensemble report](../dev/native_posterior_profile_20261110/ensemble_v1/report.json),
[locked ensemble outputs](../dev/native_posterior_profile_20261110/ensemble_v1/RNA_LOCKED.json).

Full suite:477Python passes,1GPU-only skip,8existing HiGHS-option forwarding
warnings,153.16s. All4Rust tests pass. New tests add21cases to the preceding
456 without changing any existing test. This increment uses CPU/SciPy for the
new estimator; no new GPU parity or biological coverage claim follows.
[Python report](../dev/native_posterior_profile_20261110/validation_23158642/pytest.txt),
[Rust report](../dev/native_posterior_profile_20261110/validation_23158642/cargo_test.txt).

## Accuracy evidence transported by exact prediction identity

The retained single-reference profile has exposed regression meanPearson
0.7648167082 versus strongest pseudobulk0.7302613671. The fixed eight-reference
window ensemble has0.7704112303 versus identically averaged pseudobulk0.7327616523.
These are prior48-patient/24-24split outcomes with23 informative evaluation
correlations, not new independent tests. Native equality establishes delivery
of that estimator, not a new performance experiment or restored statistical
identifiability. Prior calibration, finite support and uncertainty limitations
remain. The broader69/249-family promotion gates fail. Reference averaging
reduces loss recall.602747→.540476; this package option does not repair that
tradeoff or justify default promotion. The just-tested robust-median alternative
was rejected, not silently substituted.
[Original profile evidence](POSTERIOR_REGION_AVERAGE_2026-11-09.md),
[reference mean and failed uncertainty](REFERENCE_BAGGED_PROFILE_2026-11-10.md),
[rejected median](ROBUST_REFERENCE_PROFILE_2026-11-10.md).

Computationally:the specified port and invariants pass. Statistically:RNA
working-model interpretation only; conditional uncertainty/superiority remain
failed or unestablished. Biologically:no new independent validation. These are
separate assessments. Public exposure is optional descriptive functionality,
not a change of bulk calls, clone search/confirmation/merging or scientific criteria.

## Execution, resources and preservation

| Job | Role | Account/QOS | Actual CPU / RAM / time limit | Outcome |
|---|---|---|---|---|
|23158635|private dependency copy|innovation|1 /8GB /1h|complete0:0,4m12s|
|23158642|noneditable wheel/full suites|nslab|2 /32GB /1h|complete0:0,3m18s|
|23158657|all48central installed profiles|innovation|1 /8GB /2h|complete0:0,3m28s|
|23158658|all39eight-reference ensembles|nslab|1 /4GB /1h|complete0:0,13s|

The build initially requested4cores/2h, then reduced to2cores on the SAME
original005node after confirmed Resources pending with2free cores. The retained
per-node/task minima were synchronized. A1hbuffer remains far above prior
build duration. Every scontrol change/state is retained; no running work was
restarted/canceled and no source/test was changed for scheduling.
[Exact adjustments](../dev/native_posterior_profile_20261110/resource_adjustments.json).

Analyses use capturedruntime003 on005, a new private Conda prefix, immutable
source, neutral cwd and protected frozen paths. Dependency versions are copied
from the existing validated prefix, with no new dependency or shared install.
Clone/full-pilot jobs keep their preceding environment/snapshots untouched.
All4jobs complete0:0. Pytest process peakRSS2,131,608KB; scheduler peak1,829,020K.
Complete analysis/setup accounting remains linked rather than treating missing
short-job scheduler RSS as actual zero memory.

Wheel SHA256 `86190a0590c725ae9a243bc45dbfea40429899cb5676c92e196e1e9e329f552f`.
Native SHA256 `e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`,
identical to the preceding validated native binary.
[Installed artifact](../dev/native_posterior_profile_20261110/validation_23158642/installed_artifact.json),
[all commands](../dev/native_posterior_profile_20261110/jobs.json),
[accounting](../dev/native_posterior_profile_20261110/accounting_final.tsv),
[resources](../dev/native_posterior_profile_20261110/resource_summary_final.json),
[preservation](../dev/native_posterior_profile_20261110/preservation_final.json).

Only one new package module and one new test file are added. README adds usage.
Research delivery runner is [native_posterior_profile.py](../benchmarks/hardening/native_posterior_profile.py).
No previous package/native/test module changes; frozen scientific data/results
and every168prior ledger values stay. Documentation records the complete effort.

The independent shared-prevalence clone experiment now passes3/5full original
cases, including first historical failure108after merging/assignment;109/110
remain. This increment does not claim those failures fixed,5%whole-specimen
error control, suitable biological membership truth or SOTA completion.
