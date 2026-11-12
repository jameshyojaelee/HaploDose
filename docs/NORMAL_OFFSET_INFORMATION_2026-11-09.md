# Normal-calibration information scaling: partial improvement, no default

Information-scaled CV reduces regional error relative to reusing a
half-normal-sample shrinkage fraction on a full normal refit. It does not remove
the bias, establish normal/query transport, or repair the original clone
recovery and error-control failures. No production change is promoted.

## Locked question and measurement model

Task `dev/normal_offset_information_20261109` follows
`NORMAL_OFFSET_SHRINK_2026-11-09.md`. Its `DESIGN.md` was saved before fresh
seed-905 observations were generated. The two task endpoints are bin-average
squared error and regional-contrast squared error against the known DNA-linked
RNA response. They are different estimands; improvement in one cannot establish
the other.

The observations are independent raw binned NB counts with known nominal
exposures, baseline probabilities, depth nuisance and dispersion. Each of
16 paired simulation draws has 60 calibration normals, 240 query cells and
40 bins. Five conditions share draws: null, gain, shared regional nuisance,
mismatched normals and tumour transcriptional state. Biological specimens: zero.
These are not observed-library, gene-level, complete-workflow or biological
validation experiments.

The proposed mechanism is a sample-information mismatch: predictive CV estimates
a half-fit correction, but the full normal refit has more information. The
competing explanation is that normal prediction itself prefers attenuation of a
regional effect, even after accounting for sample information. The decisive
observation is whether information scaling resolves regional error on fresh
paired draws without discarding the bin-error endpoint or failed conditions.

## Research implementation and independent numerical gate

For normal-cell probabilities p and mean counts m, expected NB logit information
is A' diag(w) A, with row j of A equal to e_j-p and
w_j=m_j/(1+phi_j*m_j). The research implementation uses its trace divided by
J-1 as a scalar information summary. Information is evaluated at the known
baseline plus normal nuisance, without query counts or scoring truth.

For fold information I0,I1, set Ibar=(I0+I1)/2. Normal-only opposite-half
prediction selects rho, equivalent to prior precision Ibar*(1-rho)/rho.
Training and full fits use their own I/(I+precision). The full correction is
2*rho/(1+rho). Endpoints are handled algebraically without numeric infinity.
Both NB and pooled normal-profile families use the same NB heldout scoring law.
Original fold scores remain distinct from the descriptive full refit.

This is a local isotropic approximation, not a full multivariate NB posterior.
The numerical gate verifies 12 exact rational isotropic quadratic solutions
against SciPy optimization, full Fisher matrices using independent complex-step
log-mean Jacobians at two step sizes, and additivity, duplication, gene permutation
and logit-gauge invariants. Its anisotropic quadratic counterexample has objective
0.55095 for scalar shrinkage versus 0.48090909 for the exact matrix solution.
That limitation is retained, not repaired by loosening a tolerance.

Roundoff bounds depend on operation count and magnitude. Quadratic parameter
discrepancy is bounded through objective curvature. The independent selection
audit uses explicit NB log probabilities on 33/65-point meshes with local
refinement, including boundary-adjacent intervals. It checks objective agreement
within its accumulation budget, then repeats both scientific comparisons with
the independently refined parameters. All 80 comparison orderings are stable.

Default seed-904 generation replays all 1,760 prior observation/truth/challenge
arrays exactly; prior normal selection replays all 80 records and 480 normal
arrays exactly. No prior result is overwritten.

## Fresh conditional results

Each entry is the mean across 16 simulation draws, not 16 biological replicates.
The four NB arms use identical observations, support and query fitting. No
smoothing, outcome-based reference selection or query-based shrinkage is added.

| Condition | NB correction | Bin MSE | Regional MSE |
|---|---|---:|---:|
| Null | None | 0.00085123 | 0.00007363 |
| Null | Full anchor | 0.00538206 | 0.00063668 |
| Null | Fixed CV fraction | 0.00094600 | 0.00009914 |
| Null | Information-scaled CV | 0.00110347 | 0.00012170 |
| Gain | None | 0.00114455 | 0.00014330 |
| Gain | Full anchor | 0.00675388 | 0.00071023 |
| Gain | Fixed CV fraction | 0.00132081 | 0.00016795 |
| Gain | Information-scaled CV | 0.00154983 | 0.00019019 |
| Shared nuisance | None | 0.03000303 | 0.12392827 |
| Shared nuisance | Full anchor | 0.00622491 | 0.00044380 |
| Shared nuisance | Fixed CV fraction | 0.00505729 | 0.00869652 |
| Shared nuisance | Information-scaled CV | 0.00489723 | 0.00317195 |
| Mismatched normals | None | 0.00114455 | 0.00014330 |
| Mismatched normals | Full anchor | 0.03268112 | 0.12564985 |
| Mismatched normals | Fixed CV fraction | 0.02007916 | 0.07028600 |
| Mismatched normals | Information-scaled CV | 0.02481238 | 0.09150775 |
| Tumour state | None | 0.01448409 | 0.05979514 |
| Tumour state | Full anchor | 0.01804067 | 0.06314163 |
| Tumour state | Fixed CV fraction | 0.01458608 | 0.05972893 |
| Tumour state | Information-scaled CV | 0.01468311 | 0.05985708 |

The shared-nuisance regional truth is 0.24327906. NB mean response is 0.24097232
with full anchoring, 0.33157109 with fixed-fraction CV, and 0.29349015 with
information scaling. Exact empirical decomposition attributes regional squared
mean errors of 0.00000532, 0.00779548 and 0.00252115, respectively; corresponding
between-draw variation is 0.00043848, 0.00090104 and 0.00065080. Thus residual
regional error is not just sampling variation. This decomposition is descriptive;
between-draw variation includes random baselines/exposures and counts.

The pooled baseline has the same qualitative tradeoff. Its shared-nuisance
bin/regional MSE is 0.00571150/0.01038245 for fixed-fraction CV and
0.00546626/0.00407518 for information scaling, versus
0.00683318/0.00051616 for full anchoring. With mismatched normals the
information-scaled pooled bin/regional MSE is 0.02481787/0.08817456, worse than
fixed-fraction CV at 0.01989753/0.06680127. Every paired case and all eight
methods remain in the machine-readable records.

The candidate changes both fold-specific fractions during selection and the
full-refit mapping. This comparison does not isolate the mapping's contribution
from the change in selected penalty. It supports a partial procedure improvement
under shared nuisance, not a unique causal account of the prior regional bias.

Gain and tumour-state observations and fitted profiles are exactly equal in
every arm, despite different DNA labels. Normal-only selection cannot distinguish
those biological explanations. Normal matching and the RNA/DNA bridge require
independent measurements; computational correctness cannot establish them.

## Independent evidence and disposition

All paths below are under `dev/normal_offset_information_20261109/`:

- `numerics_v1/report.json`: independent numerical gate and default replay.
- `fit_v1/{references.json,reference_arrays.npz,fitted.npz,records.json,report.json}`:
  sealed normal selection, query estimates and all metrics.
- `audit_v1/report.json`: independent replay of 960,000 count entries,
  240 scalar pooled normal fits, 320 scalar selection records (128 distinct
  inputs), 640 metric records and 80 numerical ordering checks for both endpoints.
  Fisher information is recomputed using explicit per-cell matrices.
- Query-only and truth-only challenges preserve all 480 normal arrays and
  normal-selection JSON bytes. Truth-only changes preserve all 640 fitted vectors.
- `decomposition_v1/report.json`: 1,640 exact rational per-bin/regional empirical
  loss decompositions, also checking the saved regional endpoint.
- `source_integrity.json`: all 76 existing package/native/test/README/build
  files remain identical to the initial dirty source.
- `jobs.json`, `accounting_final.tsv`, and `analysis_JOB/`: exact commands,
  source/wheel/import identities, CPU backend, runtime, peak process RSS and exits.

Jobs 23103466, 23104095 and 23104250–23104254 all completed successfully.
Each used CPU1/4GB/1h on innovation or nslab, without arrays. Maximum observed
process RSS was 196,688 KiB. Fit runtime was 10.46s; independent audit 6.51s.
Slurm's zero/missing sampled RSS is not used as a memory measurement.

The research changes are in `normal_offset_information.py`,
`normal_information_numerics.py`, and optional paths in the existing panel,
shrinkage, audit and decomposition scripts. The original normal selector and
default seed remain. Installed production code and its previously validated
395-Python/4-Rust suite are unchanged; the full suite was not rerun here.

Computational status: the stated approximation and isolated conditional experiment
pass their independent checks. Identification: normal transport and the genomic
interpretation remain unsupported in the challenged conditions. Biological
validation: absent (n=0). Reject an unqualified default and do not claim the
regional problem solved. The prior original recovery result remains 2/5;
no new complete-workflow recovery, false-split, detection-power or uncertainty
coverage claim follows. Those criteria and all frozen artifacts are unchanged.
