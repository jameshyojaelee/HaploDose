"""Reference-resampling bands for a conditional RNA depth-transport prototype."""
import argparse
import importlib.util
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta
from haplodose.clones import reference_depth_slopes
from haplodose.provenance import start_run, finish_run

spec = importlib.util.spec_from_file_location('estimation', Path(__file__).with_name('depth_reference_estimation.py'))
estimation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(estimation)
CONDITIONS = ('supported', 'extrapolated', 'narrow', 'small', 'zero', 'confounded')


def shape(expected):
    value = np.log(expected)
    return value-np.median(value)


def calculate(arguments):
    panel, out, condition, rep = arguments
    key = f'{condition}__{rep}'
    with np.load(panel/'observations_and_targets.npz') as data:
        Y = data[key+'__Y']
        query = data[key+'__query_exposure']
        null = data[key+'__null']
        gain = data[key+'__gain']
    original = np.load(panel/(key+'_predictions.npz'))
    centre = shape(original['estimated'])
    rng = np.random.default_rng(np.random.SeedSequence([2026110910, rep, CONDITIONS.index(condition)]))
    samples, residuals, indices = [], [], []
    for _ in range(199):
        index = rng.integers(0, len(Y), len(Y))
        counts = Y[index]
        depth = counts.sum(axis=1).astype(float)
        gamma = reference_depth_slopes(counts, depth, [np.arange(len(Y))], 30)
        gamma -= gamma.mean()
        expected, _, _, record = estimation.predict(counts, gamma, query, 1e-12)
        samples.append(shape(expected.sum(axis=0)))
        residuals.append(record['reference_residual'])
        indices.append(index)
    samples = np.array(samples)
    difference = samples-centre
    sup = np.max(np.abs(difference), axis=1)
    regional = np.abs(difference[:, :10].mean(axis=1))
    rank = math.ceil(.95*(len(samples)+1))
    simultaneous = float(np.sort(sup)[rank-1])
    fixed = float(np.sort(regional)[rank-1])
    # Read target means only for evaluation, never for fitting or quantile choice.
    records = []
    for event, target in [('null', null), ('gain', gain)]:
        contrast = shape(target)-centre
        region = float(contrast[:10].mean())
        records.append(dict(condition=condition, rep=rep, event=event,
            whole_profile_band=simultaneous, fixed_region_band=fixed,
            region=region, region_lower=region-fixed, region_upper=region+fixed,
            positive_region=bool(region-fixed > 0), excludes_zero=bool(abs(region) > fixed),
            null_profile_uncovered=bool(event == 'null' and np.max(np.abs(contrast)) > simultaneous)))
    np.savez_compressed(out/(key+'_bootstrap.npz'), samples=samples, centre=centre,
        resampled_indices=np.array(indices), calibration_residual=np.array(residuals))
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    panel = a.task/'panel_v2'
    run = start_run(a.out, dict(script=Path(__file__), design=a.task/'UNCERTAINTY_DESIGN.md',
        estimation=Path(estimation.__file__), data=panel/'observations_and_targets.npz',
        preceding_audit=a.task/'audit_v1/run_manifest.json'), {},
        dict(scope='Reference-only bootstrap; no clone/biological calibration'), mode='depth_reference_uncertainty')
    try:
        arguments = [(panel, a.out, condition, rep) for condition in CONDITIONS for rep in range(8)]
        workers = min(4, int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = [row for group in pool.map(calculate, arguments) for row in group]
        frame = pd.DataFrame(rows)
        frame.to_csv(a.out/'intervals.tsv', sep='\t', index=False)
        summaries = []
        for condition in CONDITIONS:
            null = frame[(frame.condition == condition) & (frame.event == 'null')]
            gain = frame[(frame.condition == condition) & (frame.event == 'gain')]
            failures = int(null.null_profile_uncovered.sum())
            n = len(null)
            summaries.append(dict(condition=condition, simulation_replicates=n,
                whole_profile_null_failures=failures,
                null_failure_upper95=float(beta.ppf(.95, failures+1, n-failures)) if failures < n else 1.,
                null_failure_upper95_simultaneous6=float(beta.ppf(1-.05/6, failures+1, n-failures)) if failures < n else 1.,
                fixed_region_null_exclusions=int(null.excludes_zero.sum()),
                fixed_region_gain_positive=int(gain.positive_region.sum()),
                median_profile_band=float(null.whole_profile_band.median()),
                median_region_band=float(null.fixed_region_band.median())))
        pd.DataFrame(summaries).to_csv(a.out/'summary.tsv', sep='\t', index=False)
        (a.out/'report.json').write_text(json.dumps(dict(summaries=summaries, bootstrap_samples=199,
            quantile_order=190, biological_specimens=0, production_changed=False,
            scientific_criteria_changed=False, scope='Conditional RNA intervals, not specimen clone false-split or power'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU installed slope routine and SciPy reference-only bootstrap')


if __name__ == '__main__':
    main()
