# What remains in the observed RNA counts when population parameters are known

The scientific goal remains genomic profiles and supported tumour clones,
validated at the claimed resolution. The previous goal turn made progress by
falsifying two repairs and completing their independent numerical checks. The
goal remains active and unachieved. This continuation diagnoses whether the
exposed clone failure is attributable to fitted population parameters or the
six excluded regional genes,without another latent-rate simulation bank.

## Exact question and assumptions

Use the SAME archived seed202611099001,600 cells and actual raw counts for the
state and two-population conditions. Both use known baseline,size10 and full
region compositional response. P has exactly240 uniformly placed shifted cells,
delta=.24327906486489861,plus Gaussian state sigma=.11918191480409052.
Q has one Gaussian-state population with the same marginal mean/variance as P.
These RNA observation laws are the amplitude oracle's specified broader
competitor;Q is not the original narrower phase9 false-split calibration null.

Compare logP/logQ on original selected genes and every raw gene in the known
affected region. Fixed-k marginalisation uses the degree240 symmetric polynomial;
the iid mixing alternative is retained separately. No latent amplitudes or cell
labels enter calculation. Delta,sigma,k,baseline and full region are oracle
inputs;this is not a usable inference method. Estimand:relative likelihood of
these fixed RNA population laws for this observed specimen. No population error
or power rate is estimated. Regional transcription can reproduce P;the score
does not identify DNA dosage or genomic clones. Biological n=0,base seeds=1.

Keep24/48/96 Gaussian quadrature and original scientific criteria. Independently
check exhaustive unbalanced fixed-count and zero-effect fixtures,long-double
recurrences,a different600-factor polynomial convolution,complete-PMF adaptive
integration at consequential cells,and raw count/baseline/gene/cell maps.
Numerical allowances derive from recurrence operation counts and convergence.
They are computational diagnostics,not confidence intervals or p-values.
[Frozen full design](../dev/count_population_oracle_20261111/DESIGN.md).

## Completed result

| Condition | Genes | Fixed-count logP/logQ | Numerical diagnostic | IID-mixture score |
|---|---:|---:|---:|---:|
| State | Selected94 | −22.459774 | 0.001258 | −0.097891 |
| State | Full100 | −22.493397 | 0.001292 | −0.096584 |
| True contrast | Selected94 | −0.601988 | 0.001985 | −0.698832 |
| True contrast | Full100 | 0.076639 | 0.002566 | −0.038028 |

The selected event retains94/100 affected genes. Knowing the parameters does
not produce a large positive likelihood contrast in this actual true specimen.
Restoring all six genes changes its sign but leaves likelihoods close. This
falsifies the sufficient-repair claim for parameter knowledge alone on the
selected measurements. It does not establish a finite-count power limit or
prove that no valid procedure could reject this observation.

[Inference] Missing measurements can contribute,but this result does not justify
promoting a larger selected region as the general repair. Fitted-model cost,
count sampling and the broader biological competitor remain distinct. A small
observed log-ratio in one draw cannot be substituted for a distribution of
oracle scores or for population detection power. Do not excuse the original
recovery/error failures or weaken their criteria. The existing noiseless-rate
power bound is not an RNA-count guarantee.

Every independent polynomial and full-PMF check passes. All12 full-size
polynomials,36 adaptive component integrals,raw counts,known baselines and
gene/cell maps agree within unchanged numerical diagnostics. No generating
truth labels/amplitudes were opened. No new null bank or controller promotion.
[Authoritative independent report](../dev/count_population_oracle_20261111/reader_23214996/results/REPORT.json).

## Executed code and resources

Jobs23214995/23214996 complete0:0 in25s/11s on nslab/nslab CPU Intel-v4,
1CPU8GB1h each,reader afterok oracle. Existing isolated Conda environment,
neutral directories,source hashes before/after,time-v resources. No arrays,
production/default changes,new dependencies,shared installs or frozen edits.
The earlier running reader23213185 and unrelated jobs remain untouched.

[Study](../dev/count_population_oracle_20261111/study.py),
[independent reader](../dev/count_population_oracle_20261111/reader.py),
[source snapshot](../dev/count_population_oracle_20261111/SNAPSHOT.json),
[exact submission commands](../dev/count_population_oracle_20261111/submissions.json),
[oracle resources](../dev/count_population_oracle_20261111/run_23214995/resources.txt),
[reader resources](../dev/count_population_oracle_20261111/reader_23214996/resources.txt).
