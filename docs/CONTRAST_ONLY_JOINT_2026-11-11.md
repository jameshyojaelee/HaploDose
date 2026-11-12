# Contrast sharing with unrestricted expression state

The preceding coupled model fixes clean low-depth whole-profile loss but fails
the full gate, passing 28/32 RNA and 55/64 simultaneous tests, with high-depth state
regressions. This matched ablation retains adjacent discrete contrasts while
freeing full-resolution continuous expression-state loadings. It tests that
specific cause without a parameter sweep or weaker criteria.

[Frozen scientific contract](../dev/contrast_only_joint_20261111/DESIGN.md),
[native constrained implementation](../dev/contrast_only_joint_20261111/restriction.py),
[complete regression runner](../dev/contrast_only_joint_20261111/study.py),
[independent reader](../dev/contrast_only_joint_20261111/reader.py).

Every count, fold, reference, pairing offset, selection rule and original
baseline remains. Reuse unchanged no-factor fits exactly and require bit-exact
replay of all their predictions and held-out scores. Factor-enabled fits use
the old coupled fits as feasible nested parents, with all fixed native starts
and existing quadrature/optimizer settings. The working model's state-zero
response depends on Gaussian independent state; a better RNA fit cannot identify
DNA dosage, genomic membership or biological confidence. Biological n=0.

Qualification23219504 reaches the report writer after the numerical assertions,
then fails on `KeyError: objective`: reused coupled no-factor fit dictionaries
store the training objective in their nested coarse metadata. Their native
parameters are valid. Source_v1, snapshot, log and resources remain unchanged.
Source_v2 changes one report expression to evaluate the objective directly with
the same native model. No numerical assertion or scientific criterion changes.
The full batch remains gated on an actually completed qualification.

[Failed qualification](../dev/contrast_only_joint_20261111/numerical_23219504/output.txt),
[first immutable source](../dev/contrast_only_joint_20261111/SNAPSHOT_v1.json).

Qualification23219524 completes0:0 in0.60s with59924KiB process RSS. All eight
independent full-PMF/reduced-gradient fixtures pass, along with half-step,
label/factor symmetry, identity embedding and feasible-nested-fit checks.
The largest gradient discrepancy is recorded in the numerical report;budgets
derive from its operation bounds. Both failed and corrected qualification
sources remain. All338 prior ledger values and91 protected identities match.

Eight full blocks23219535–23219542 and dependent reader23219543 are submitted.
Each block requests1CPU8GB3h;reader1CPU8GB2h,all account/QOS nslab,CPU partition,
without arrays. Each computes two offsets of every fixed case:290 cases/580
phase bundles. No scientific result follows until all locks and independent
checks pass. Original32 RNA and64 simultaneous comparisons remain unchanged.

[Passed qualification](../dev/contrast_only_joint_20261111/numerical_23219524/results/REPORT.json),
[current source identities](../dev/contrast_only_joint_20261111/SNAPSHOT.json),
[all batch commands and dependencies](../dev/contrast_only_joint_20261111/submissions.json),
[exact case catalogue](../dev/contrast_only_joint_20261111/case_jobs.json).

The first eight blocks23219535–23219542 all stop1:0 after5–6s at the first
fold-model writer: reused JSON dictionaries already contain `key`, so adding
that field again raises `TypeError`. No complete phase lock or scored result
exists from that batch. Reader23219543 cannot satisfy its failed dependencies
and is cancelled while pending. Source_v2 and every failed output remain.

Source_v3 removes only that redundant metadata key when loading old model
dictionaries. It changes no parameter, likelihood or selection. A new bounded
preflight23219587 uses1CPU4GB1h,nslab CPU,testing the first lexical seed at both
depths of all three original conditions,each at both offsets. It executes the
entire count-only writer and independent phase reader without calling truth
targets or scoring scientific performance. The full batch must await this pass.
All original numerical/scientific checks remain;failed jobs are not hidden.

[Failed batch example](../dev/contrast_only_joint_20261111/run_0_23219535/output.txt),
[second immutable source](../dev/contrast_only_joint_20261111/SNAPSHOT_v2.json),
[preflight command](../dev/contrast_only_joint_20261111/preflight_submission.json),
[preflight implementation](../dev/contrast_only_joint_20261111/preflight.py).

Preflight23219587 completes0:0 in2m08s,109588KiB peak process RSS and73624KiB
scheduler-sampled RSS. All12 writer/independent-reader cases pass, including
both depths and both offsets of homogeneous, clean two-population and state
conditions. No truth targets or scientific outcomes were opened. The unchanged
no-factor fold scores and posterior arrays are exact. Source_v3 is frozen.

Full blocks23219647–23219654 and dependent reader23219655 replace the failed
adapter batch,with identical resources and scientific work. The old catalogue
and its checksum are preserved as `case_jobs_failed_v2.*`;old submissions retain
their original job paths. Only metadata adapters changed. No inference failure
has been hidden or scientific criterion relaxed. Full results remain pending.

[Integrated preflight report](../dev/contrast_only_joint_20261111/preflight_23219587/results/REPORT.json),
[successful-source batch commands](../dev/contrast_only_joint_20261111/submissions_v3.json),
[failed-batch catalogue](../dev/contrast_only_joint_20261111/case_jobs_failed_v2.json).

## Completed result: reject promotion

All290 cases/580 phase bundles lock and all eight successful-source fitting jobs
complete0:0. Reader23219655 completes0:0 in4m32.09s,116000KiB process RSS and
95416KiB scheduler RSS. It verifies10440 exact score sums,2320 held-out model
pairs,3480 posterior arrays and3441 adaptive full-PMF probes. Every unchanged
no-factor score and point is exact. Numerical validation passes.

The scientific gate fails:28/32 original RNA and52/64 simultaneous both-target
comparisons pass. This is below the original joint model's30/32+57/64 and the
coupled model's28/32+55/64. Freeing expression state is not a sufficient repair.

| Whole-profile MSE,48 fixed seed units/condition | Original joint | Coupled | Contrast-only | Pseudobulk |
|---|---:|---:|---:|---:|
| Clean/300,RNA |.00431896|.00335300|.00429474|.00413259|
| Clean/300,state-free |.00407038|.00384453|.00385499|.00413259|
| Clean/3000,RNA |.000340596|.000617286|.000270567|.00348458|
| Clean/3000,state-free |.000167622|.000612131|.000919615|.00348458|
| State/300,RNA |.00725041|.00767485|.00705693|.01619652|
| State/300,state-free |.00429165|.00393385|.00409065|.00440117|
| State/3000,RNA |.000563710|.00242837|.000897691|.01554463|
| State/3000,state-free |.000366723|.000762995|.000754100|.00375644|

The high-depth RNA recovery relative to coupled fits supports a cost from
restricting state loadings. But clean high-depth RNA improvement accompanies
worse state-free profiles, and clean low-depth RNA loses the pseudobulk gain.
These outcomes reject interpreting improved RNA prediction as sufficient
dosage-profile improvement. The contrast-sharing restriction may route genomic
signal into the continuous factor; a saved-fit decomposition is required before
claiming that mechanism. No additional threshold or penalty search follows.

All12 failed simultaneous tests remain:

| Target | Condition/depth | Endpoint | Comparator | Resolved wins/48 |
|---|---|---|---|---:|
| RNA |Clean/300|Whole|Pseudobulk|23|
| RNA |Clean/300|Whole|Global|35|
| RNA |State/300|Minority|Beta expected logit|26|
| RNA |State/3000|Whole|Beta expected logit|23|
| RNA |State/3000|Minority|Individual|18|
| State-free |Clean/300|Minority|Global|22|
| State-free |Clean/300|Minority|Beta expected logit|9|
| State-free |Clean/3000|Minority|Individual|1|
| State-free |Clean/3000|Minority|Beta expected logit|3|
| State-free |State/300|Minority|Pseudobulk|8|
| State-free |State/300|Minority|Global|7|
| State-free |State/300|Minority|Beta expected logit|11|

The35/48 global comparison passes the original32 family but fails the jointly
protected64 family;both original rules remain. Homogeneous RNA error remains
above pseudobulk. No biological clone, uncertainty, default or SOTA claim changes.
All21 jobs in this experiment are terminal:11 complete0:0,9 failed1:0 retained,
and one impossible dependent reader cancelled while pending. Successful fits
take9m46–14m23,peak process RSS89228–96424KiB. All use account/QOS nslab.

[Complete independent result](../dev/contrast_only_joint_20261111/reader_23219655/results/REPORT.json),
[all21 job outcomes](../dev/contrast_only_joint_20261111/accounting_complete.tsv),
[reader resources](../dev/contrast_only_joint_20261111/reader_23219655/resources.txt).
