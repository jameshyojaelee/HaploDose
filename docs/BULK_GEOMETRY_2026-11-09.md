# Genomic geometry does not repair the exposed bulk-loss failure

Matched within-block distance scaling leaves mean comparison loss recall
unchanged at0.500183, slightly increases neutral-window FPR and slightly reduces
correlation. Actual spacing does not outperform the fixed shuffled-spacing
control. Retain this negative result; do not promote a new default or perform a
more elaborate distance-model search on these exposed outcomes.

## Scientific contract and decisive comparison

The estimand is regional natural-log relative RNA effect under the fitted pooled
Poisson-lognormal working model, with latent four-state RNA dosage. Reference,
programme shifts, expected counts, beta, centre and variance remain fixed to
exact reproduced HaploDose fits. This is a conditional transition-geometry
ablation, not a nuisance-refitted complete alternative workflow.

The hypothesis is that smoothing by bin number suppresses loss evidence where
measurements are physically sparse. The main competing computational explanation
for an apparent improvement is weaker overall smoothing, rather than useful
genomic geometry. Biological RNA regulation/reference mismatch remains a
competing explanation for DNA-shaped RNA profiles; this experiment cannot rule
it out. Independent DNA at suitable resolution remains necessary for genomic
interpretation. No biological replication or independent evaluation is added.

Use all48 exposed raw-STAR tumour patients and20 normal reference donors, the
original24/24 calibration/comparison split, original genes and callable support.
DNA is the same uncertain GRCh37 ABSOLUTE reference, gene anchored to hg38 RNA
at approximately10Mb windows. This is not interval liftover or new DNA truth.
No RNA predictions, priors or nuisance fits use DNA values. Prediction for every
RNA arm completes before DNA is loaded for scoring.

## Geometry controls, fixed before execution

Each bin is located at the median start of its retained genes. These measured
loci are monotone; long overlapping genes can make interval midpoints unsuitable.
Within each declared chromosome/centromere/gap block, successive distances are
divided by their block mean. Total process time therefore equals the original
number of edges. This controls total smoothing exposure rather than quietly
making long blocks less regularized. It does not assert a universal biological
breakpoint rate per base pair.

The ideal equal-exit transition matrix uses the existing state prior and
p_change. Its principal matrix logarithm is checked to be a real Markov generator
with nonnegative off-diagonal rates; no clipping or ridge repairs invalid
transitions. Exponentiating that generator by scaled distances preserves the
block endpoint transition matrix by the semigroup identity. Zero distances have
identity transitions. Existing reset boundaries remain fixed.

Compare the exact legacy kernel, an ideal constant-transition numerical control,
actual distance scaling and one fixed within-block distance permutation
(seed20261109). The permutation preserves each block's weight distribution and
process time but removes association with actual spacing. It is a mechanism
control, not a permutation test or sufficient evidence for a causal mechanism.
The original comparator family and strong RNA-only naive partition baseline
remain included, under identical truth/callable support.

## Independent computational checks

Preflight23067836 verifies reversible spectral powers independently of matrix
logarithm/exponential construction; maximum discrepancy8.88e-15 is below the
operation-scaled binary64 allowance8.53e-14. Integer matrix powers and semigroup
identities agree. Exhaustive256-state-path enumeration includes zero distance
and a chain reset; its independent float64 recursion differs by2.22e-16, and
the JAX recursion by4.15e-8, within the existing1e-4 HMM criterion. These are
numerical checks, not biological thresholds or formal floating-point enclosures.

Panel23068287 reconstructs counts and all48 tumour fits exactly. Legacy per-bin
posteriors/effects, original truth/support and patient loss metrics reproduce
exactly. Every geometry arm is also checked against independent float64
block-wise forward/backward calculations using the same emission inputs.

Audit23068763 independently reconstructs all1131 bin anchors from sorted gene
starts and all39 blocks from chromosome/block identities. Spectral matrix
products check matched block endpoints, including shuffled distances, within
operation-scaled arithmetic allowances. All1440 patient recall/FPR/precision
checks agree; exact Fraction calibration ranks agree. Recalibration of the
float64 reference predictions changes no callable-window decisions in any arm,
despite small floating-point threshold differences. No tolerance is relaxed.

## Descriptive results and paired populations

Each score uses the original patient-equal neutral-window5% calibration rule on
the24 calibration patients, with conservative strict-threshold ties. That target
does not establish a5% achieved comparison FPR or a whole-specimen clone error
guarantee. The latter is a different estimand and remains unvalidated.

| Exposed comparison role | Mean loss recall | Mean neutral FPR | Mean Pearson |
|---|---:|---:|---:|
| Legacy posterior |0.500183|0.057043|0.720291|
| Ideal constant control |0.500183|0.057043|0.720291|
| Actual distance |0.500183|0.057312|0.720238|
| Shuffled distance |0.511172|0.056389|0.720702|
| Strong naive RNA-only partition |0.399679|0.041244|0.730261|

There are24 comparison patients; loss recall is defined for13, neutral FPR and
correlation for23. Patients, not windows, are the replication units. Shared
reference donors remain a dependence source. No new p-values, biological
confidence intervals, independent superiority or calibration claims are made.

Actual versus legacy spacing has paired recall difference0, FPR difference
+0.000269888 and correlation difference-0.000052825. Precision difference is
+0.000021649 over the same14 patients with defined precision in both arms.
Actual versus shuffled spacing has lower recall, higher FPR and lower
correlation, with a small positive precision difference. One permutation does
not support a general inferential comparison, but supplies no positive evidence
for the proposed geometry explanation in this bounded experiment.

Naive partitioning has higher correlation and lower FPR but lower loss recall;
neither procedure dominates across tasks. Precision is defined for14 patients
in each method, but those are different sets. Their paired precision comparison
uses only the same13 patients, giving distance-minus-naive difference0.006084.
Do not compare unpaired macro precisions as though they shared a population.
Undefined metrics, support counts and all patient rows remain explicit.

All RNA effect arms are in natural-log units. The preserved original naive
scores use log2(CPM+1) reference contrasts; the strong partition baseline uses
their natural-log-scaled version. Calibration is score-specific; raw amplitudes
from different definitions are not compared as interchangeable DNA effects.

## Decision, limits and reproducibility

The ideal constant control excludes a numerical explanation for these point
results. Matched process time excludes weaker total prior exposure as an
unreported intervention. Actual geometry supplies no demonstrated recovery
benefit here. Do not promote this variant, tune a physical rate to the existing
DNA, or extrapolate the negative result to every possible distance-dependent
model. Nuisance refitting and a globally uniform per-bp process were not tested;
neither is automatically justified by a more elaborate mathematical model.

Implementations are `benchmarks/hardening/genomic_geometry.py`,
`bulk_geometry_panel.py` and `bulk_geometry_audit.py`. Designs, immutable source
snapshots, exact nuisance replays, emissions, anchors/distances/permutation,
window predictions, independent references, paired effects, hashes, environment
and actual CPU backend records are under `dev/bulk_geometry_20261109`.

| Job | Account/QOS | Stage | Seconds | Process peak RSS KiB | Exit |
|---|---|---|---:|---:|---:|
|23067836|innovation|preflight|15|471564|0|
|23068287|nslab|48-patient conditional panel|45|974948|0|
|23068763|innovation|independent audit|3|85608|0|

Production code, original tests and frozen artifacts are unchanged. The full
package suite was not rerun for these diagnostic-only additions; all new
executable checks and analyses ran on compute nodes. Original recovery2/5,
clone false-split/power failures, uncertainty failures and missing biological
null panel remain unresolved. No SOTA or full-goal completion claim follows.
