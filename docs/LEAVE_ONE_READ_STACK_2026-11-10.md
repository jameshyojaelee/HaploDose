# Test whether half-depth selection explains the remaining RNA profile loss

The near-full-depth correction does not establish a sufficient accuracy fix.
It improves all four locked state cases, improves two of four clean CNA cases
and worsens all four homogeneous cases. These are 12 exposed mechanism cases
from six base seeds at paired depths, not a population-performance result.
No full comparison gate or package default is promoted.

## Mechanism and independent evidence

The preceding stack learns weights from N/2-read forecasts and uses those
weights for N-read target forecasts. The new path scores every possible
removed read at N−1 reads, weighted by the count of that read category.
For fixed multinomial probabilities, the exact identity is
`P_N(y) * y_j/N = p_j * P_(N-1)(y-e_j)`.
Independent rational enumeration verifies all 18 supported three-category
cases at total three. This supports a fixed forecast's reserved-read risk;
it does not make an adaptively selected risk unbiased.

The estimand is conditional regional RNA probability and its exported
relative RNA log odds. The existing six experts, priors, chromosome-respecting
partitions, 96 EM iterations, cell folds, normalization and within-group
allocation are unchanged. Source counts are the existing fixed-total
multinomial simulation. Its conditional splitting identity does not validate
biological NB sampling or identify genomic dosage from RNA-equivalent state.

Weighted SLSQP is followed by independent pair-coordinate minimization.
The convex gradient gap must reach its matrix-specific arithmetic estimate;
otherwise the fit abstains to the pool. A known asymmetric optimum is checked
with a parameter allowance derived from its directional gradient uncertainty
and curvature bound 4/9. Canonical row ordering, exact target-fold exclusion,
an unchanged peer forecast and positive normalized compositions pass. Two
independent rational derivative fixtures differ by zero and 1.11×10⁻¹⁶,
within their 1.07×10⁻¹⁴ short-operation allowance.

## Every locked mechanism case

Values are whole-profile squared RNA log-odds errors. Each row retains the
same generated data and scoring population as the saved refined half-depth
stack. All 96 corresponding baseline score records reproduce exactly.

| Condition | Seed | Reads | Half-depth selection | Leave one read |
|---|---:|---:|---:|---:|
| Homogeneous | 9641000 | 300 | 0.00163778 | 0.00167728 |
| Homogeneous | 9641001 | 300 | 0.00115114 | 0.00128242 |
| Homogeneous | 9641000 | 3,000 | 0.00012297 | 0.00012470 |
| Homogeneous | 9641001 | 3,000 | 0.00014525 | 0.00014875 |
| Two populations | 9642000 | 300 | 0.00415739 | 0.00428017 |
| Two populations | 9642001 | 300 | 0.00406418 | 0.00400822 |
| Two populations | 9642000 | 3,000 | 0.00059146 | 0.00059042 |
| Two populations | 9642001 | 3,000 | 0.00065422 | 0.00066286 |
| Two populations plus state | 9643000 | 300 | 0.00901471 | 0.00872917 |
| Two populations plus state | 9643001 | 300 | 0.00855991 | 0.00823426 |
| Two populations plus state | 9643000 | 3,000 | 0.00231661 | 0.00220493 |
| Two populations plus state | 9643001 | 3,000 | 0.00202745 | 0.00202032 |

All 1,080 regional fits satisfy the numerical criterion. The mixed clean
results do not support depth mismatch as a sufficient explanation of the
remaining accuracy failures. They also do not prove that its expected effect
is zero in every regime. Weak information, model misspecification, expert
nonidentifiability and normalization remain alternatives. A full performance
run is not used to turn this pilot into a successful gate. All eight baselines,
minority/state-free errors, contrasts and fit records remain available.

## Execution and limits

Job 23189745 completes `0:0` in 1m02s on nslab CPU with one CPU, 8 GB
and one hour requested. Its timed study uses 49.95 seconds and 87,552 KiB
peak RSS. Installed wheel/native/backend identities and all source hashes
are checked. The independent derivative reader completes as 23190127;
its first serialization failure is retained in the linked transfer report.
Biological n=0. Original recovery assertions, 5% specimen criteria,
uncertainty requirements, frozen outputs and production files are unchanged.

[Design](../dev/leave_one_read_stack_20261110/DESIGN.md),
[inference code](../dev/leave_one_read_stack_20261110/study.py),
[weighted numerical solver](../dev/leave_one_read_stack_20261110/weighted.py),
[submission](../dev/leave_one_read_stack_20261110/submission.json),
[complete report](../dev/leave_one_read_stack_20261110/run_23189745/results/report.json),
[all scores](../dev/leave_one_read_stack_20261110/run_23189745/results/scores.tsv),
[all fits](../dev/leave_one_read_stack_20261110/run_23189745/results/decisions.json),
[independent derivatives and paired recount](../dev/scone_regional_shrink_20261110/summary_23190127/results/SUMMARY.json),
[real-data transfer and attempt history](SCONE_REGIONAL_SHRINK_2026-11-10.md).
