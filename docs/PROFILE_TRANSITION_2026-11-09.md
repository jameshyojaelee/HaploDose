# Measure the applied profile update, not an unused returned normal dose

The prior conditional-objective diagnostic scored the raw `shared_profiles`
return. Refinement actually fixes the normal-component dose to zero and can
discard components before applying the profile update. The corrected diagnostic
records that exact transition before changing offset, dispersion, depth slopes
or programme scores. On unchanged production inference, applied profile updates
decrease conditional NB Q in 31/67 updates, versus 22/67 for the raw return.
Original recovery remains 2/5. No production change follows from this alone.

## Locked contract and independent checks

Task `dev/profile_transition_20261109` saves `DESIGN.md` before rerunning the
five original immutable recovery inputs. All six original criteria per seed
are retained. Biological specimens: zero.

The estimand is conditional cell-level NB log-likelihood change, holding
responsibilities, nuisance and dispersion fixed. When a component is dropped,
old and new doses use exactly the same retained responsibility columns, without
renormalizing their mass. This avoids comparing different latent populations.
It does not measure dropping components' contribution to the complete observed
mixture likelihood.

The NB RNA model differs from the pooled PLN observation model used for
segmentation. Lower NB Q can reflect segmentation regularization, shared-dose
constraints or differing working likelihoods. It alone establishes neither a
production defect nor a biological failure; improving it does not identify DNA.

Research instrumentation in `offset_nb_objective.py` adds default-off
`trace_applied`. It retains original raw-profile and offset records and observes
the applied transition at `_update_offset` entry, after normal-row and keep
decisions. Contexts save pre-filter responsibility/dose arrays, retained indices,
old/new doses and native responsibility precision.
`offset_objective_panel.py --trace-applied` saves both traces.

Independent `profile_transition_audit.py` requires:

- Exact replay of all 134 historical raw-profile/offset rows after accounting
  for inserted record indices, and all 30 original recovery criteria.
- One applied transition per raw transition, with matched model/component
  identities, fixed zero normal dose and explicitly aligned retained indices.
- Exact equality of retained arrays to saved pre-filter subsets. Responsibility
  roundoff uses native precision and component count; float64 conversion does
  not improve their original precision.
- Independent full SciPy NB logpmfs for 45 worst contexts selected by minimum
  Q within model/kind without scoring truth, plus bounded first/last-cell
  Decimal80 likelihood-ratio checks. Budgets use operation count and magnitude.

All checks pass. The diagnostic reproduces prior inference while correcting
which proposed profile transition the likelihood claim describes.

## Observations and disposition

| Quantity | Updates | Decreases beyond numerical budget |
|---|---:|---:|
| Raw returned profiles | 67 | 22 |
| Applied retained profiles, normal dose fixed zero | 67 | 31 |
| Subsequent production offset | 67 | 33 |

These are dependent refinement updates, not independent specimens or failure
probabilities. Minimum raw/applied profile gains are -35.49665430 and
-35.37749667. This comparison does not establish that the unused normal row
explains every difference; component filtering is also recorded.

Seeds 20261031 and 20261107 pass. Seeds 20261108/109/110 fail exactly the same
criteria as before. The last two still return one S1 tumour group and fail the
original assignment and chr12 difference criteria. This measurement correction
does not recover a missing clone.

Historical 26/69 and 26/63 profile-decrease counts for conditional-ratio/NB-offset
research arms also describe raw returns. Their applied counts were not run here.
Prior recovery failures and offset calculations remain unchanged. Withdraw an
applied-step interpretation of those raw profile counts.

Under `dev/profile_transition_20261109/`, `source_v1/` is the execution snapshot,
`legacy_v1/` contains all five fitted outputs, traces and contexts, and
`audit_v1/report.json` contains independent checks. `source_integrity.json`
verifies all 76 existing package/native/test/README/build files remain identical
to initial dirty source. Frozen results and original tests remain unchanged.

Jobs 23105338 and 23105339 complete on compute nodes with the readonly validated
wheel/environment: innovation CPU4/8GB/2h, then nslab CPU1/8GB/1h, no arrays.
Actual process peaks are 1,524,812 KiB and 228,144 KiB; runtimes 80.95s and
17.05s. `jobs.json`, `accounting_final.tsv` and `analysis_JOB/` preserve commands,
backend/import identities, resources and exits. Static AST checks pass. The
full production suite was not rerun because production code/tests are unchanged.

Computational status: executed-transition recording and independent replay pass.
Identification and biological validity: this diagnostic supplies neither.
Retain the corrected measurement, without replacing segmentation merely to raise
NB Q or weakening recovery, 5% false-split/power or uncertainty criteria. Those
scientific failures remain unresolved.
