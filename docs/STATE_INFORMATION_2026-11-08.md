# A continuous-state alternative limits clone detection even before count noise

Under one explicitly specified competing RNA-state model, the best possible
level-5% specimen test has detection power at most0.1118242, with at least95%
Monte Carlo confidence conditional on the stated laws and checked numerical
enclosures. This is a600-cell oracle calculation, not a universal biological
limit, a HaploDose performance estimate, or proof that the three original
thinning assertions are impossible. Those failures and criteria remain intact.

## The assumption being tested

The prior simulation fixed a single-population state distribution. This diagnostic
asks what happens when its regional mean and variance are not independently
identified. Keep600 cells, RNA contrast delta=0.6*log(1.5), and within-class
Gaussian state SD0.25. Compare:

| Law | Unobserved biological interpretation | Observed centred regional amplitude |
|---|---|---|
| P | Two CNA populations, exactly300 cells each, uniformly permuted labels | Means minus/plus delta/2, independent Gaussian state SD0.25 within class |
| Q | One CNA population with a continuous transcriptional state | Gaussian mean0, variance0.25²+delta²/4 |

Both laws have the same per-cell mean and variance. Q represents a stronger
specific nuisance alternative than the earlier fixed-mean/fixed-variance null;
it is not substituted into the frozen experiment. Unlike an independent mixture
approximation, P retains the original fixed balanced cell counts. Arbitrary cell
order carries no class information. The parameters were fixed before outcomes.

The oracle observes noiseless amplitudes and knows both distributions. Under the
earlier expression-only generator, the common depth law and NB counts depend on
class through the total regional amplitude. Counts therefore cannot make this
P-versus-Q distinction more informative than observing those amplitudes directly.
This argument does not cover additional independent DNA/allele measurements or
an externally identified state distribution that rules out Q.

## Decisive calculation

The Neyman–Pearson likelihood ratio is the optimal specimen test for these two
specified laws. The balanced P density sums over all possible balanced label
assignments, using an elementary-symmetric-polynomial dynamic program in log
space. It is not a fitted two-Gaussian clustering or a classifier trained on
known labels. Q has a direct Gaussian density.

There are20000 independent simulated null specimens and20000 alternative specimens;
biological n is zero. A lower97.5% confidence bound on Q's true95th likelihood-ratio
quantile uses its18939th ordered score. The exact failure probability is.0237817;
the next rank would give.0256235 and fails the prespecified.025 allocation.

The lower cutoff intentionally provides an upper envelope for optimal level-.05
power. It is not proposed as a usable5% decision rule. At that cutoff,2149/20000
alternative specimens are possible rejections. A one-sided97.5% exact binomial
upper bound, rounded outward, is0.11182413156542002. The two failure allocations
give at least95% joint confidence. Numerical score bounds are carried in the
conservative direction; no alternative sample is ambiguous at the boundary.

Computational falsifiers are explicit: exhaustive densities, exchangeability,
zero-effect identity, independent precision or binomial identities disagree.
These checks pass. The finding instead limits statistical identification and
power under the specified competing law. It does not establish that real tumours
follow either law, or that an RNA state truly explains a particular DNA-supported
clone. Same-cell DNA or independently justified state measurements would address
that biological question.

## Interpretation for development

[Inference] Removing membership loss can help the locked recovery cases while
remaining insufficient to distinguish genomic groups from plausible state
variation. In this particular regime, merely reducing sequencing noise cannot
produce high clone-detection power under a5% constraint against Q. More cells,
larger effects, independent state constraints or genomic evidence change the
problem; this result does not supply a universal detection threshold.

Do not tune a more permissive merge or score threshold to force recovery, promote
the failed joint procedure, or declare all RNA subdivisions false. Retain
candidate-only interpretation and abstention where genomic distinctions lack
identifying evidence. The5% biological guarantee still requires independent
specimens and suitable truth;40000 oracle simulations do not supply that panel.

## Reproducible evidence

Jobs22997738/22997739 complete in1m34/1m36 on innovation/nslab CPU nodes. Main
source is candidate_v1 under `dev/state_information_20261108`. All latent vectors,
labels, scores, parameters, seeds and numerical allowances are retained.
Summary22997802 is preserved; final summary/audit22997969 takes25s and adds
outward rounding and independent density verification without resimulation.

Small-n exhaustive SciPy enumeration, zero-effect and cell-permutation identities,
and600-cell extended-precision checks pass. A separate posterior-count polynomial
convolution verifies six actual extreme/threshold-adjacent vectors. Its maximum
discrepancy is2.85e-13, inside the recorded operation-count allowances, at most
3.36e-9 across the simulations. These are checked floating-point error estimates,
not formal interval arithmetic. Independent exact-integer polynomial sums verify
the calibration rank and bracket the power root;60/80-digit calculations agree.
The final report is `dev/state_information_20261108/summary_v2/report.json`.
