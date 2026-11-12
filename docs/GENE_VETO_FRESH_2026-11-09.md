# Fresh gene-count workflows reject the clone false-split claim

The frozen aligned-pool candidate recovers the planted two-profile chr12
contrast in98/100 new simulated specimens, but reports additional profiles in
41/100 one-CNA specimens with continuous expression state. Its simultaneous
false-split upper bound is53.45%, above the unchanged5% target. Reject promotion.
The comparator fails similarly. The separate weight-unit failure independently
blocks the candidate. Production remains unchanged; exposed5/5 recovery does
not override these negative results.

## Target, observation model and alternatives

This is a fresh conditional simulation experiment, not biological evaluation.
Each specimen has600 query cells from one or two generating CNA populations,
600 reference cells across three types/splits, and4401 source genes. Common
chr3 gain/chr8 loss, depth3000, gene phi0.1, beta0.6 and reference/depth/mismatch
generation follow the existing simulator. References and source nominal query
means are paired across conditions. Query means are redrawn as independent
gene gamma-Poisson counts after prescribed regional changes and normalization
to the original nominal exposure. These are synthetic raw integer counts.

The three locked conditions are a matched one-CNA null, a one-CNA state null,
and two CNA populations with the same state variation. State is a continuous
Gaussian log-rate effect over the first half of chr12. The two-population
contrast adds delta=.6*log(1.5) to a randomly chosen40% of query cells. State SD
is delta*sqrt(.6*.4), matching the between-population log-effect SD of that
prespecified60:40 contrast. It is a dimensioned challenge, not a clinical cutoff.
The two-clone condition contains both state and the CNA contrast; this panel
does not variance-match its total variation with a stronger state-only null,
nor cover a binary state observationally equivalent to a CNA.

The primary estimand is the whole-specimen probability of reporting an
unsupported additional non-reference profile within the one established CNA
population. The separate power estimand is probability of exactly two final
non-reference profiles and only chr12 differences in the implanted contrast.
This does not demand exact breakpoints/effect sizes and is not cell-membership
accuracy. The event is read from public clone tables after terminal validation;
it is not a selected-parent score or a post hoc count of well-assigned groups.

A competing explanation is discretization of continuous expression state.
It is distinguished here by known generating CNA labels, not by calling RNA
groups genomic truth. Suitable corresponding DNA with documented sensitivity
would be needed to distinguish it biologically. Both procedures still know
true gene size10. Nuisance, full-count profiles, membership weights and adaptive
search remain fitted; conditional thinning does not make that inference independent.

## Locked procedures, units and replication

Both procedures run the complete pipeline from feature/reference selection
through nuisance fitting, initialization, all proposals, confirmation, gene
vetoes, EM, merging and stopping. Both use coherent gene-native count splits
and the existing membership-weighted confirmation/merge stages. Only the
candidate aligns the gene-veto pools with actual regional weights. Original
model seed, support rules, nominal sign cutoff, pseudocount and all numerical
thresholds stay fixed. No new calibration bank or tuned operating point is used.

`dev/gene_veto_fresh_20261109/DESIGN.md` predates new outcomes. A disjoint
six-specimen pilot(2/condition,12 paired fits) checks execution and scoring and
is excluded from primary evaluation. The primary panel contains100 independent
simulated specimens per condition, paired across arms and conditions, totaling
600 fitted workflows. There are100 independent generating seed units, not600
biological replicates; biological n=0. No draws are added after inspecting bounds.

The four null upper endpoints(two laws x two procedures) and two recovery lower
endpoints form one fixed Bonferroni family6, simultaneous95% one-sided exact
Clopper-Pearson bounds, per-endpoint alpha=.05/6. These are a new experiment's
bounds; prior experiment families and failed bounds remain unchanged. A null
upper bound must be<=5%. Power is assessed simultaneously without inventing
a universal minimum. Missing/failed fits cannot disappear from denominators;
there are no execution failures here.

## Results, including negative evidence

| Condition | Procedure | Primary events /100 | Simultaneous bound | Interpretation |
|---|---|---:|---:|---|
|Matched one-CNA null|Comparator|0|upper4.6747%|passes this conditional null|
|Matched one-CNA null|Aligned candidate|0|upper4.6747%|passes this conditional null|
|Continuous-state one-CNA null|Comparator|39|upper51.4404%|fails5% criterion|
|Continuous-state one-CNA null|Aligned candidate|41|upper53.4487%|fails5% criterion|
|Two CNA populations + state|Comparator|97|lower90.0760%|two-profile/chr12 recovery only|
|Two CNA populations + state|Aligned candidate|98|lower91.6455%|two-profile/chr12 recovery only|

The candidate adds one recovered contrast and two state-null false distinctions
relative to its comparator. No significance or superiority claim is made from
those paired changes. Its state-null outputs have1/2/4 non-reference profiles in
59/14/27 specimens; the comparator has61/13/26. These are final reported profile
counts, not a substitute definition based on MAP assignments. A future diagnosis
of sparsely assigned profiles must retain this frozen event definition and result.

Mean assignment accuracy over all600 query cells, penalizing reference/unresolved
assignments, is99.885% for both under the matched null,95.6633%/95.4917% for
comparator/candidate under the state null, and72.7917%/72.8000% under two clones.
Mean tumour-call fractions are99.885%,99.8383%/99.8333%, and99.8833%/99.8833%.
These are descriptive specimen means, not calibrated membership probabilities
or patient-level accuracy estimates. High contrast recovery does not establish
accurate cell assignments. The constant one-group baseline has zero null splits
and zero exact two-profile recovery by definition; it cannot qualify as clone
recovery despite suppressing errors.

The matched-law success and state-law failure localize a scope failure of the
observation/search procedure. They do not alone identify which stage creates
each false group, prove that only model misspecification matters, or establish
that every possible RNA-based procedure fails. Earlier state-null rejections
are retained alongside this distinct gene-count experiment.

## Independent checks and primary measurement chain

The research hook gains only an explicit caller-provided case-ID allowlist;
the original five IDs remain its default. Original replay23120316 reproduces
all five candidate recovery outcomes, exact cell/difference tables, all30 criteria
and all15 split arrays. No prior scientific output changes. Pilot23120317 and
audits23120375/23120414 pass. The latter explicitly regenerates source references
and nominal query means before primary jobs begin.

Primary chunks23120418–23120429 all complete. Audit23120430 independently checks
all300 specimen-condition input chains and600 workflow scores. It regenerates
all100 source reference-count and nominal-mean identities, reconstructs state/label
streams, and checks every generated query count via the specified gamma-Poisson
streams. Independent log-space normalization checks mean construction. Truth-free
metadata, source feature/cell order, raw integer CSR, region mapping and unchanged
reference observations agree across conditions. Both arms' actual gene split
arrays match exactly for every paired input. Each fit records actual CPU backend.

Public clone/difference/cell tables independently reconstruct every primary event,
assignment, call fraction and group count. Exhaustive cell-label matching checks
the producer's Hungarian matching. Each exact-beta bound agrees with independent
binomial-tail root inversion. Transitive input/output/source/wheel/native hashes
are retained in manifests. This is computational correspondence, not a biological
truth assay or independent interval-coverage validation.

## Execution, preservation and next decisions

All17 CPU jobs complete0:0 on innovation/nslab, no arrays or interactive QOS.
Original/pilot requests4CPU/8GB/2h; two pilot audits1CPU/8GB/1h; twelve primary
chunks4CPU/8GB/3h; final audit1CPU/8GB/2h. Final audit scheduler runtime6m38s,
process6m35.75s. Maximum process RSS across task3379156KiB, below3.23GiB.
Exact commands, execution inventory, backend, dependency status and per-job
resources are in jobs.json, accounting_final.tsv and resource_summary.json.

Source/code are gene_veto_fresh.py, gene_veto_fresh_audit.py and the explicit
case-ID extension in gene_thinning_hook.py under benchmarks/hardening. Immutable
source_v1/audit_source_v1/final_source_v1 and execution records are in the task.
All76 protected package/native/test/README/build files remain unchanged. Full
package tests are not rerun for research-only code; prior isolated wheel evidence
is395 Python/4 Rust passes,1GPU-only skip. Frozen inputs, old results and all
historical decisions remain intact.

Do not promote either tested procedure as5%-controlled clone inference. Preserve
the strong conditional recovery and both scientific failures. The separate
[unit falsification](GENE_VETO_UNITS_2026-11-09.md) requires a versioned measurement-
unit repair, not a favorable rescaling or reanalysis of these outcomes. Next trace
representative state-null multi-profile outputs through initialization, support,
assignments and terminal merge, including sparsely assigned profiles, before
changing the procedure. Any future correction gets new frozen data and the same
error/power requirements; it cannot overwrite this test or relabel these100
specimens as an unopened evaluation set.

Gene-dispersion estimation, membership/nuisance/selection uncertainty, repeated-
workflow coverage, profile/bulk/spatial accuracy, suitable biological null and
same-cell DNA evidence, and fixed-competitor superiority remain unresolved.
Neither computational success nor more simulations supplies missing biological
replication or makes the full manuscript/SOTA objective achieved.
