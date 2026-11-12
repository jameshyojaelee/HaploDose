# Ordered pooling in the joint state/profile model

The scientific goal remains accurate genomic dosage and supported tumour clones
across RNA modalities, established against strong comparators and suitable
independent measurements. The cached-model diagnostic rejects simply forcing
two populations: low-depth minority error falls while whole-profile error grows,
and oracle memberships leave substantial profile error. This experiment tests
the smallest ordered regularization of the same joint count likelihood.

[Completed cause trace](JOINT_MEMBERSHIP_TRACE_2026-11-11.md),
[frozen design](../dev/adjacent_joint_profile_20261111/DESIGN.md),
[execution source identities](../dev/adjacent_joint_profile_20261111/SNAPSHOT.json).

## Statistical contract

Aggregate adjacent pairs of regions within each supplied genomic block. Fit the
unchanged conditional-multinomial mixture/Gaussian-factor model to those counts,
then expand to the original regions using MLE within-pair fractions from the
training RNA population. The full likelihood factors into the pooled model
and conditional within-pair compositions. Thus the restriction makes population
contrasts and continuous loadings equal within each pair while preserving a
full-resolution shared RNA trunk. It reduces fitted contrast/loading dimensions
without a fitted penalty or chosen copy-number amplitude.

Fit both pairing offsets and average point estimates with fixed weights 1/2.
Offsets reset at every block boundary. No DNA outcome chooses a grid or weight.
The toy panel has 24 ordered regions in one declared block; pooling two bins
does not establish a physical genomic resolution. The two restricted models
can disagree at genuine boundaries. Within-pair state or genotype variation
violates their shared-composition assumption. Matched controls and both targets
must reveal that cost.

Embed the pooled parameters into the original full-feature count model before
computing expected log odds. This preserves the exact observation likelihood;
it does not invert a mean probability. Keep RNA and state-free estimands
separate. The latter sets the Gaussian coordinate to zero and remains a model
counterfactual. Zero-mean Gaussian state independent of population and adequate
shared loadings are identifying restrictions. RNA-equivalent transcriptional
populations remain competing explanations for a CNA label; suitable independent
DNA must distinguish them. Biological n is zero.

Phase-specific responsibilities are diagnostic RNA mixtures. No component labels
are averaged into a combined genomic membership output. Public support stays
30, and the 24-cell minority cannot be called a supported genomic clone. This
multinomial experiment does not replace the heterogeneous gene-NB false-split
panel or establish uncertainty calibration, absolute copy number or ancestry.

## Fixed full test and controls

All 288 original exposed cases retain their raw counts, references, folds,
24-region evaluation support and 48 generating seed units per scenario paired
across depth. Add the two already fixed state-only controls. Both offsets yield
580 count-only bundles; all lock before truth scoring. No new seed, feature,
depth, optimiser, initialisation, quadrature or hyperparameter sweep.

Keep every original 32 RNA test and simultaneous 64 both-target test, including
all four original baselines and exact one-sided win-probability bounds. Full
qualification still requires every criterion. Boundary-overlap decisions do not
become wins. Passing permits separately frozen transfer, not a default or broad
scientific validation claim.

Matched pseudobulk and individual controls use the same pair sums, conditional
RNA fractions, original-region expansion and equal-offset average. Their
half-count binomial-margin interpretation stays explicit. Descriptive all-cell
fits share all-cell fractions; predictive folds use training cells only. Retain
the original joint/no-factor profiles and all historical baselines. A benefit
explained by the shared operator does not support a special joint-model mechanism.

## Compute and independent qualification

Numerical job **23218617** completed `0:0` on `ne1dc6-003`. Its eight combinations
of offset, component count and factor status pass count conservation, block
isolation, likelihood factorization, full PMFs, responsibilities and both point
estimands. Scalar checks verify the matched controls. The identity partition
reproduces parameters and points exactly. Core optimiser/model source remains
byte-identical to its independently qualified version.

[Numerical report](../dev/adjacent_joint_profile_20261111/numerical_23218617/results/REPORT.json),
[qualification source](../dev/adjacent_joint_profile_20261111/qualify.py),
[pooling and embedding](../dev/adjacent_joint_profile_20261111/pooling.py).

Eight nslab-account/nslab-QOS CPU AMD-v6 jobs **23218626–23218633** request one
CPU, 8 GB and three hours each. Blocks contain 36 original cases, with both
state controls added to the last block. Independent reader **23218634** requests
one CPU, 8 GB and two hours after all eight succeed. No arrays, installations,
shared-environment changes or running-job modifications. All eight fitting jobs
completed `0:0`. Reader 23218634 failed an integration-resolution assertion at
original case 144 after verifying cases 0–143. No partial scientific result is
accepted. Its failed log and original source remain unchanged.

[Exact submissions](../dev/adjacent_joint_profile_20261111/submissions.json),
[case catalogue](../dev/adjacent_joint_profile_20261111/case_jobs.json),
[independent reader](../dev/adjacent_joint_profile_20261111/reader.py).

The reader reconstructs training fractions and embedded parameters independently,
checks exact predictive sums, frozen component choices, full PMFs, high-order
posterior points, adaptive integrals, both controls and exact-tail decisions.
Numerical allowances include quadrature/refit convergence and phase errors
conservatively; they are not statistical confidence intervals. Source hashes,
immutable snapshots, actual node/environment, runtime and memory remain recorded.
The 335 prior ledger values, 91 protected identities, frozen results and all
scientific criteria stay unchanged. The full objective remains unachieved.

## Reference-integrator failure and separately versioned correction

The old adaptive check divides by the unrestricted saturated full-profile PMF.
The pooled model can have a very small conditional within-pair likelihood that
is constant over the integrated state. Its scaled integrand then falls below
the absolute-error resolution, even though the model density is well defined.
The assertion correctly prevents accepting that unresolved check. This is a
reference-computation problem, not a result about model performance.

The independent numerical correction keeps the scientific reader byte-identical
and every fitted parameter, point prediction, mask, endpoint, tolerance and
64-test decision unchanged. For each discrete component, log density `g` has
`g''(t) = -1 - N Var_p(loading) <= -1`. At any approximate mode `t0`, its global
upper envelope is `g(t0) + g'(t0)^2/2 - (t - t0 - g'(t0))^2/2`. Remove that
constant scale, integrate around the derived centre, and bound tails by the
same declared Gaussian-tail tolerance. This cancels tiny additive likelihood
constants. Mixture probabilities and errors are accumulated in log space.

Qualification **23218968** requests one nslab CPU, 4 GB and one hour. It must
reproduce the actual original failure, pass all four small component/factor
cases against old/full-PMF references, check the score with independent complex-
step derivatives and repeat the original consequential-cell checks. Only a
passing report permits the full reader. No fitting rerun, science-criterion
change, skipped cell, epsilon increase or favorable-seed substitution.

[Failed reader output](../dev/adjacent_joint_profile_20261111/reader_23218634/output.txt),
[failed reader resources](../dev/adjacent_joint_profile_20261111/reader_23218634/resources.txt),
[exact correction design and derivation](../dev/adjacent_joint_reader_20261111/DESIGN.md),
[new source identities](../dev/adjacent_joint_reader_20261111/SNAPSHOT.json),
[qualification submission](../dev/adjacent_joint_reader_20261111/numerical_submission.json).

## Completed result: reject general promotion

Qualification23218968 passes, including the reproduced original failure,
independent complex-step score checks and all consequential-cell probes.
The separately versioned full reader23218992 then completes0:0 in4m52.47s,
107076KiB peak process RSS; scheduler sampling reports85620KiB. No fits or
predictions were rerun. Every original tolerance and scientific criterion remains.

The model passes28/32 original RNA tests and55/64 simultaneous both-target tests.
The complete gate fails. With48 fixed seed units per heterogeneous stratum,
clean low-depth whole RNA MSE improves from original joint0.00431896 to0.00335300,
below pseudobulk0.00413259 and matched-coarsened pseudobulk0.00411539. This fixes
the specific low-depth whole-profile failure. But high-depth state-plus-profile
RNA MSE rises from0.000563710 to0.00242837, and high-depth state-free MSE from
0.000366723 to0.000762995. The original all32 high-depth simultaneous successes
are not preserved. Homogeneous RNA error improves, but still exceeds pseudobulk.

All nine failed simultaneous criteria remain:

| Target | Condition/depth | Endpoint | Comparator | Resolved wins/48 |
|---|---|---|---|---:|
| RNA | State/300 | Minority | Beta expected logit |29|
| RNA | State/3000 | Whole | Beta expected logit |0|
| RNA | State/3000 | Minority | Individual |9|
| RNA | State/3000 | Minority | Beta expected logit |29|
| State-free | Clean/300 | Minority | Global |21|
| State-free | Clean/300 | Minority | Beta expected logit |7|
| State-free | State/300 | Minority | Pseudobulk |16|
| State-free | State/300 | Minority | Global |13|
| State-free | State/300 | Minority | Beta expected logit |19|

The independent reader verifies10440 exact score sums,2320 held-out model pairs,
3480 posterior arrays and3444 adaptive full-PMF probes. Both fixed state-only
controls selectK1 at both offsets. No combined clone assignments or biological
error estimates are emitted. All14 jobs in this trace/ordered-model chain are
terminal:13 complete0:0, with the initial failed reader1:0 retained. All use
account/QOS nslab. Neither additional compute nor selected endpoint gains
overrides the failed gate.

The next matched ablation frees only continuous state loadings while keeping
discrete contrasts shared within pairs. This tests whether coupling state to
genomic adjacency caused the high-depth regression; that explanation is not yet
established. No new penalty or threshold search follows.

[Completed independent report](../dev/adjacent_joint_reader_20261111/reader_23218992/results/REPORT.json),
[passed numerical qualification](../dev/adjacent_joint_reader_20261111/numerical_23218968/results/REPORT.json),
[all14 job outcomes](../dev/adjacent_joint_reader_20261111/accounting_complete.tsv),
[contrast-only ablation](CONTRAST_ONLY_JOINT_2026-11-11.md).
