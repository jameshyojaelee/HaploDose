# Exact patient weighting repairs a loss-threshold boundary defect

The development scorer's threshold could change when identical neutral windows
were duplicated within each calibration patient. This changes the representation
without changing the intended probability distribution. Exact rational weighting
repairs the defect. The existing real-data point thresholds are unchanged;
the completed 20,000-draw nested replay also preserves point effects and loss
probabilities. Million-draw loss precision remains under investigation.

## Estimand and decisive counterexample

The calibration distribution gives each patient with usable neutral windows equal
total mass, spread uniformly over that patient's windows. Repeated bootstrap
patients contribute repeated mass. The threshold is the smallest observed score
whose cumulative mass reaches the declared 5% quantile. Calls use the strict rule
`score < threshold`, preserving conservative treatment of ties.

Twenty patients with constant scores 0 through 19 define the same distribution
whether each score appears once or in many identical windows. The exact lower
5% quantile is 0. The old calculation returned 1 for 47 of 100 window-replication
counts. Fractional weights accumulated in floating point crossed an exact mass
boundary inconsistently. This establishes a computational invariance defect;
it does not establish that the old empirical false-positive fraction exceeded
5% in this fixture. Neither calculation establishes population error control.

The correction in `benchmarks/hardening/profile.py` uses an integer common
denominator for patient weights and exact cumulative-mass comparisons. The
declared decimal probability is interpreted as a rational number, so 0.05 is
1/20. Bisection over the observed score support finds the same lower quantile
without jitter or a numerical tolerance. Patients without usable neutral windows
remain excluded from that calibration distribution as before; an entirely empty
calibration set remains unresolved.

## Validation and limits

Compute jobs 22999297 and 22999336 retain the before/after counterexample:
47 failures become zero across the same 100 representations. Job 22999380 passes
all 30 focused scorer tests. New checks include 600 heterogeneous distributions,
independently enumerated with `Fraction`, repeated patient indices, unequal
coverage, missing observations and reordered patients/windows. Original test
assertions remain unchanged.

Job 22999399 checks all 225 thresholds across nine frozen reference fits and
25 methods. Independent exact cumulative probabilities bracket every threshold.
No point threshold changes; corresponding patient-level scores match exactly.
The replay retains 21,600 historical/corrected patient-score rows. No new model
fit, patient, truth construction or comparator is introduced.

Bootstrap weights differ from the point weights. The corrected nested replay
22999447 therefore retains the original seed, patient/reference resampling and
153-comparison family in new artifacts. The original factorial analysis and its
precision jobs continue against their immutable source. Do not treat the point
replay as proof of bootstrap parity, or promote loss results before the corrected
replay and numerical precision checks finish.

The corrected nested job22999447 completes in38m36s, peak497836KiB. All point
effects and loss probabilities reproduce the original factorial analysis; its
historical loss-probability change list is empty. Independent prefix audit22999741
checks20000 exact resampling identities,51 loss contrasts,204 scalar full-scorer
effects and188 interval endpoints. Collector check23000148 also passes. These
checks admit the fixed remaining980000 draws, now running as25 independent CPU
jobs23000360–23000384. Collector23000385 and independent audit23000411 follow.
This validates the correction's computational consistency on this panel; final
precision and biological interpretation remain separate gates.

The inference remains conditional on exposed development patients, eight normal
reference resamples, the stated DNA truth and common callable support. Exact
quantile arithmetic supplies neither suitable biological replication nor verified
DNA copy number, cell membership, or confidence-interval coverage.

## Evidence

Artifacts and exact executed sources are under
`dev/loss_quantile_identity_20261108`: `probe_v1`, `corrected_v1`, `tests_v1` and
`replay_v1`. The corrected nested analysis is separately specified in
`dev/gene_factorial_20261108/EXACT_THRESHOLD_REPLAY.md`. Frozen historical outputs
are retained; this is a scorer correction, not a changed scientific target.
