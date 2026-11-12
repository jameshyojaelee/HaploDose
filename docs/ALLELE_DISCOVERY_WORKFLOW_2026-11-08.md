# Withheld SNPs stay out of the complete conditional discovery fit

A research workflow now physically withholds reserved SNP observations before
the entire specimen fit, rather than only masking them at the final proposal.
Both existing cell observation models pass exact discovery-output invariance
under confirmation-only perturbations. Reserved summaries and a separate
descriptive all-SNP fit leave the frozen discovery outputs unchanged. This
repairs the demonstrated computational dependency in the conditional binned
workflow; it is not a production default or a calibrated clone method.

## Scientific contract and executed intervention

The estimand is computational dependence of the full discovery result on
designated reserved allele counts, at fixed RNA observations, reference,
dispersion and SNP metadata. The working likelihood is the existing NB RNA
model plus phased allele likelihood, with either legacy binomial cell evidence
or explicitly supplied beta-binomial concentration 2. These are deterministic
allele fixtures and one fixed NB RNA draw, not biological replication or
samples from a population performance study.

The prior counterexample held memberships fixed. This experiment calls the
actual complete `fit_specimen` search, including initialization, allele-based
group retention, EM, nuisance updates, proposals, agglomeration, split/merge,
allele vetoes and terminal validation. The research intervention constructs a
fresh allele object using only twenty metadata-selected proposal SNPs. Thus all
of those discovery consumers see only proposal allele counts. It starts with
the fixed reference concentration 50, without inheriting concentration or phase
parameters already fitted to the complete query. Genomic SNP/bin correspondence
and distances remain explicit.

The original minimum observed SNP, depth, cell, posterior and segment thresholds
are unchanged. The original internal even/odd logic operates within discovery
data, so fewer observations remain available. Computational separation is not
assumed to preserve power. No gene-veto input or spatial graph is represented;
this is a conditional binned workflow, not the full public raw-input pipeline.

## Exact invariance and nontrivial partitions

Each cell model has three odd-SNP perturbation scenarios. Proposal alternate
counts are fixed at 10, 12 or 16 of 20 reads in three forty-cell groups. Reserved
counts are 10, 19 or 1 of 20. Both models also receive two absent-proposal-coverage
cases and a repeated baseline after unrelated fits. All RNA inputs are identical.

| Cell model and workflow | Odd balanced | Odd same orientation | Odd reverse orientation |
|---|---:|---:|---:|
| Binomial, current all-SNP path |40/80|0/120|0/120|
| Binomial, discovery-only path |40/80|40/80|40/80|
| Concentration 2, current all-SNP path |2/118|0/120|0/120|
| Concentration 2, discovery-only path |65/55|65/55|65/55|

These are MAP group sizes, including the reference-like group. They are not
counts of validated tumour clones. The discovery-only cases have two groups
above the unchanged 30-cell minimum, so the invariance check is not simply
an always-single-group result. With absent proposal coverage, both models
return MAP sizes 1/119 and retain only one group above that minimum. This is
reported rather than called successful recovery.

The discovery-only path preserves every recorded fitted array, nuisance
parameter, search log, selected interval, candidate definition and SNP-index
mapping exactly. Both repeated baselines are exact. Current all-SNP discovery
changes under reserved perturbations, including its label-invariant final
partition. No floating tolerance is used for these decisions.

## Reserved evidence and separate descriptive refitting

After discovery is saved, the workflow reads reserved SNPs only to sum counts
for frozen MAP groups and discovery-selected regions. It returns raw count
summaries, not significance, an accepted split, a new group or a confidence
interval. Missing eligible pairs or regions retain an explicit empty-output
reason. No reserved-count threshold is chosen here.

A fresh all-SNP allele object then fits descriptive states to the frozen
memberships. It cannot overwrite discovery parameters or the reserved summaries.
File and array hashes are checked after both stages. The independent audit
recomputes seven reserved contrasts with scalar cell/site sums and verifies
descriptive pooled counts with independently implemented dense multiplication.
The latter allowance scales with the 120-cell floating-point accumulation;
it is not a scientific decision threshold.

Role isolation does not imply statistical independence: shared cell state,
correlated ASE, phase, reference uncertainty and selection still matter. The
strongest biological alternative remains transcriptional/allelic cell state
rather than a DNA difference. Corresponding independent DNA and suitable
replication would distinguish those explanations; byte-identical discovery
artifacts cannot. No coverage, power, false-split or clone-identity claim is
established by these fixtures.

## Verification and decision

Jobs 23002991 and 23002992 complete in 23s and 21s, using innovation and nslab
CPU allocations respectively. Process peak RSS is 692940 and 691236 KiB.
Independent audit 23003015 completes in 2s with 78988 KiB process RSS. It verifies
18 full discovery fits, reconstructs 304 fitted arrays, checks label-invariant
partitions and all saved fingerprints, and independently sums seven reserved
contrasts. Both observation models, absent proposal evidence, fitting order and
separate descriptive outputs pass the declared computational checks.

Production source, existing tests, defaults and scientific criteria are unchanged.
The wheel remains the version validated by 348 Python and four Rust tests with
one GPU-only skip; that unchanged suite was not rerun. The original thinning
result remains 2/5 and the three failures remain unresolved. The joint error/power
gate remains failed; biological null n remains zero. The unsupported joint
allele/expression-null configuration remains rejected.

Retain the workflow as a tested research implementation. Before integration,
audit public preprocessing and stable SNP roles, reference/bias filtering and
donor/phase provenance. Any scientific promotion additionally requires the
same observed/null channels and complete search, with a frozen fresh evaluation
of false splits and power. Do not substitute these fixtures for that experiment.

Implementation: `benchmarks/hardening/allele_discovery_workflow.py` and
`allele_discovery_workflow_audit.py`. Immutable design, source, all input matrices,
discovery artifacts, reserved/descriptive outputs, hashes and actual backend/
resource records are under
`dev/allele_discovery_workflow_20261108/{binomial_v1,beta_v1,audit_v1}`.
