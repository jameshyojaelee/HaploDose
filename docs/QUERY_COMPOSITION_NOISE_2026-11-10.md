# Test gene-derived bin noise in membership, not only event proposals

**Completed; correction rejected.** Job23184729 exits0:0 in3m42s. All fifteen
fits complete; independent numerical checks, all five exact sham controls and
source/proposal gene-fold equality pass. Recovery remains4/5, so the unchanged
5/5continuation gate fails. No larger null/power validation or public integration
is warranted for this unchanged correction. Production remains2/5; its three
hard failures and the biological false-split objective remain unresolved.

| Seed | Baseline / candidate complete recovery | Baseline assignment | Candidate assignment |
|---|---|---:|---:|
|20261031|pass / pass|.836032389|.838383838|
|20261107|pass / pass|.738095238|.702970297|
|20261108|pass / pass|.756916996|.756916996|
|20261109|fail / fail|.602000000|.600397614|
|20261110|pass / pass|.814000000|.809236948|

Assignment uses the original tumour-call-conditional denominator; it is reported
with all other original criteria, not substituted for exact two-group recovery.
No previously passing binary criterion is lost and no S3group is added, but109
still fails exactly-two-groups, assignment>.65 and chr12-only differences. In107,
the lower assignment accuracy is a real descriptive regression even though it
stays above the original threshold. All failed and passing measurements remain.

The actual109chr12proposal reaches confirmation in both arms. Mean membership
margin changes.658452138→.644478418, and both have zero initial weighted-HMM
segments: failure still precedes merging. S1 initial native median dispersion
is.023508781 versus proposed.024491880. These are fitted working-model noise
estimates, not independently known technical variance. Updating bin noise by this
formula does not supply the missing initial confirmation evidence. Do not lower
posterior support or add draws to make this intervention pass.

The candidate exercises15initial models and68noise-update calls. All fifteen
initializations save complete fitted fractions/noise/gene maps/rows and original
noise; the candidate keeps its values fixed through all68calls. The sham calls
the original update. Scheduler peak1345440KiB; measured process peak1604476KiB,
wall3m34.92s. The full installed package/GPU suite was not rerun for research-only
code; the newly executed numerical, exact-control and recovery checks are the
relevant verification. Biological specimens=0, no false-split guarantee.

[Final report](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/probe_v1/report.json),
[independent numerical checks](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/probe_v1/numerics.json),
[exact sham controls](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/probe_v1/SHAM_LOCKED.json),
[complete manifest](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/probe_v1/run_manifest.json).
Earlier preparation checkpoint below retains its original context.

Job23184729 submitted on CPU005,innovation account/QOS,1CPU8GB90h. It is pending
at the initial checkpoint; no recovery improvement or calibration result yet.
The experiment tests all five original recovery cases under baseline, exact sham
and candidate, after independent numerical checks. The continuation gate is5/5,
no lost passing individual criterion and no added S3group. No array, new package
dependency, installation change, public default or frozen-result modification.

The existing gene-aware event proposal experiment already tested both exact gene
likelihoods and moment-matched bin proposals. Both remained4/5 and failed109.
Repeating those unchanged would not address the remaining cause. This experiment
instead sets the noise used by the full membership/search path from discovery
gene composition and reference-derived gene dispersion, holding that noise fixed
instead of using hard-membership residuals and median shrinkage. The same older
research workflow is the baseline, not a newly validated method. Production still
recovers2/5 original cases; the research baseline previously recovered4/5 but failed
state-null control. No biological claim follows from a possible recovery rescue.

For independent gene NB counts, the bin variance is M+sum(m_g^2/r_g). With common
within-bin fractions p_g=m_g/M, effective bin dispersion is sum(p_g^2/r_g).
The discovery-only fraction estimate can differ from reference composition while
remaining invariant to a common bin-level rate shift. This is a moment-matched
working NB law, not exact closure. Correlated gene state, cell-dependent fractions
and nontransportable reference noise can invalidate it. Other existing nuisance
fits use full counts, so this experiment does not establish independent pipeline
confirmation. Independently corresponding DNA and biological replication remain
necessary for clone identity. Biological n=0.

Independent checks use exact Fraction calculations, scalar gene-variance sums,
SciPy NB moments, gene permutation, regional scaling, thinned variance and explicit
failure on unsupported bins. Floating allowances derive from operation count and
magnitude. The confirmation-array helper check addresses only that helper's input
contract; it does not certify workflow leakage isolation. Exact sham arrays,
tables and logs must pass before the candidate executes. Original source/proposal
gene folds must be equal across arms; downstream gene-veto row selections can
change with the fitted populations and are not equated to fixed observations.

Immutable source_v1 preserves all953files of the previous corrected-production
research snapshot and adds this driver/design. The launcher verifies installed
Python/native/wheel identities and all source hashes, uses the existing isolated
environment read-only and executes from a neutral compute directory. That
environment's prior498Python/4Rust passes are historical evidence, not a new full
suite run. New executable numerical and recovery checks are all on compute.
The initial metadata-only preparation used the unavailable command `python`;
it failed before any source copy or submission and was retried with `python3`.

[Frozen scientific design](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/DESIGN.md),
[research implementation](/gpfs/commons/home/jameslee/Haplodose/benchmarks/hardening/query_composition_noise.py),
[exact submission](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/jobs.json),
[source preservation](/gpfs/commons/home/jameslee/Haplodose/dev/query_composition_noise_20261110/snapshot.json).
