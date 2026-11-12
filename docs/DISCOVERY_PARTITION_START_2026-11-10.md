# Test a discovery-only initializer before blaming merge support

COMPLETE and rejected:23172453 exits0:0 in1m38s. Exact fixture enumeration,
permutation/translation/scale and rate anchoring pass. All three old discovery
objectives, nine predictive arrays and three guard scores reproduce EXACTLY;
confirmation count/row/gene reconstruction and conservation pass. No new start
qualifies for selection at any of the three orders. All six new local fits
report success. Successful optimization is not successful recovery.

| Order | Saved objective | New positive-separation local fit | Mirrored start |
|---|---:|---:|---|
|24|−.150172706944|−.150079139117|Same saved optimum within arithmetic uncertainty|
|48|−.150172719712|−.150079144569|Same saved optimum within arithmetic uncertainty|
|96|−.150172719712|−.150079144569|Same saved optimum within arithmetic uncertainty|

Objectives are minimized. The different local fit is worse at every order;
selecting it just to change the population interpretation is unjustified. The
saved−.06970696647confirmation score remains resolved negative under the same
.00107005305numerical diagnostic. This rejects this particular initializer as
a sufficient correction, not every possible optimizer or observation model.
It does not certify a global optimum or prove nonidentifiability. No seed,
threshold, initializer sweep, recovery expansion or new power panel follows.

Measured process peak179576KB;the scheduler reports103052KB. Both measurements
are retained with their scopes. All919executed source files remain exact and
all193older ledger values remain unchanged in the194-field completed checkpoint.
A no-op documentation patch missed its context after checkpoint writing and
made no edit;the report/resource addition then applied normally. No execution,
scientific result or original output was rerun or altered for that edit failure.

[All fits and checks](../dev/discovery_partition_start_20261110/probe_v1/report.json),
[discovery-before-confirmation lock](../dev/discovery_partition_start_20261110/probe_v1/DISCOVERY_MODELS_LOCKED.json),
[completed manifest](../dev/discovery_partition_start_20261110/probe_v1/run_manifest.json),
[execution resources](../dev/discovery_partition_start_20261110/analysis_23172453/resources.txt).

Earlier preparation checkpoint follows with its original context.

The phase8replicate2 saved discovery model assigns about2.66of600expected
cells to its smaller component and predicts confirmation worse than the null
(score−.069707). This precedes the later35MAP/25.066effective-cell merge loss.
Proposal-prior profiling already failed to restore support. This bounded probe
tests whether a missed mixture optimum explains the upstream failure.

The observation model remains the same gene-count NB mixture with a shared
Gaussian log-rate state, depth-specific centres and a common profiled mixture
proportion. The estimand is a distinguishable regional RNA population mixture.
Continuous state remains the strongest biological alternative; independent
corresponding DNA is needed to identify genomic clones. No new specimens or
independent evaluation are added.

Implement one new deterministic initializer plus its label mirror. Solve the
one-dimensional two-mean cut on discovery-only depth-centred log amplitudes;
use its separation and residual variance only as optimization starts. Preserve
the null's mean rate analytically. The same exact NB/state objective, null,
L-BFGS-B settings, quadrature orders24/48/96 and confirmation rule remain.
No truth labels or known simulator parameters enter fitting.

Independent Fraction enumeration over all64binary assignments verifies the
initializer fixture, alongside permutation/translation/scale and rate-anchor
identities. Recompute the saved discovery objectives EXACTLY. Keep the old fit
unless a successful new fit improves beyond the existing operation/magnitude
budget. Lock every discovery selection before reconstructing confirmation
counts from the archived full-minus-proposal CSR matrices. Require exact row,
gene, discovery and count-conservation identities, all9old predictive arrays
and all3old scores. A better discovery objective alone cannot qualify.

The existing resolved-positive confirmation rule is the scientific gate. A
failure rejects this initializer as a sufficient correction without a new seed,
threshold, initializer sweep or calibration expansion. A pass permits separately
frozen recovery/control testing, not adoption or an error-rate guarantee.

| Job | Account/QOS | Headers | Submission outcome |
|---|---|---|---|
|23172453|nslab|CPU4CPU8GB2h,original005/runtime003|accepted; result pending|

All917prior source files remain exact in the919-file snapshot. Static AST
parsing passes; numerical and scientific checks have not yet completed. The
old validated private wheel is read only. No package, native, dependency,
default, frozen result or original criterion changes. Existing clone jobs
23167849/23169730/23170760continue unchanged.

[Frozen design](../dev/discovery_partition_start_20261110/DESIGN.md),
[research implementation](../benchmarks/hardening/discovery_partition_start.py),
[source preservation](../dev/discovery_partition_start_20261110/snapshot_v1.json),
[exact command](../dev/discovery_partition_start_20261110/jobs.json),
[upstream evidence](PHASE8_REMAINING_TRACES_2026-11-10.md),
[rejected prior correction](PROPOSAL_PREVALENCE_MERGE_2026-11-10.md).
