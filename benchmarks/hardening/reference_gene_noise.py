"""Full recovery with gene splitting driven by reference-only noise estimates."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from haplodose import pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('reference_noise_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def estimator(records, depth_strata=1):
    if depth_strata not in (1, 5):
        raise ValueError('Use the original estimator or prespecified depth quintiles')
    cache = {}
    def size(info):
        frame = sys._getframe(1)
        while frame is not None and frame.f_code is not pipeline._run_single.__code__:
            frame = frame.f_back
        assert frame is not None
        loc = frame.f_locals
        used, meta = loc['used'], loc['meta']
        key = (str(loc['counts']), tuple(used.column))
        if key in cache:
            return cache[key]
        rows = np.flatnonzero(loc['fit'])
        assert len(rows) >= 20 and meta.iloc[rows].role.eq('reference').all()
        assert meta.iloc[rows].ref_split.eq('fit').all()
        X, gbin = pipeline.gene_loader(loc['counts'], used, rows)(np.arange(len(rows)))
        X = X.toarray().astype(float)
        assert np.array_equal(gbin, info['gbin'])
        L = np.asarray(loc['Lbin'][rows], float)
        assert np.array_equal(X.sum(axis=1), L)
        groups = (meta.dataset+'|'+meta.ref_group).to_numpy()[rows]
        if depth_strata > 1:
            stratified = np.empty(len(groups), dtype=object)
            for group in sorted(set(groups)):
                selected = np.flatnonzero(groups == group)
                edges = np.quantile(np.log(L[selected]), np.linspace(0, 1, depth_strata+1)[1:-1])
                strata = np.searchsorted(edges, np.log(L[selected]))
                stratified[selected] = [str(group)+'|depth'+str(s) for s in strata]
            groups = stratified
        numerator, denominator = 0., 0.
        details = []
        for group in sorted(set(groups)):
            keep = groups == group
            y, exposure = X[keep], L[keep]
            assert len(y) >= 2
            pair_exposure = exposure.sum()**2-np.dot(exposure, exposure)
            assert pair_exposure > 0
            rate_square = (y.sum(axis=0)**2-(y*y).sum(axis=0))/pair_exposure
            assert np.all(rate_square >= 0)
            expected_square = float(rate_square.sum()*np.dot(exposure, exposure))
            factorial = float(np.sum(y*(y-1)))
            numerator += factorial-expected_square
            denominator += expected_square
            details.append(dict(group=str(group), reference_cells=int(keep.sum()),
                factorial_moment=factorial, squared_mean_moment=expected_square))
        phi = numerator/denominator
        if not np.isfinite(phi) or phi <= 0:
            raise ValueError('Reference data do not identify a positive common dispersion; no oracle fallback')
        value = 1/phi
        cache[key] = value
        records.append(dict(counts=str(loc['counts']), training_rows=rows.tolist(),
            training_cell_ids=meta.cell_id.to_numpy()[rows].tolist(), genes=X.shape[1],
            estimated_phi=float(phi), estimated_size=float(value), groups=details,
            depth_strata=depth_strata,
            query_rows_used=0, holdout_reference_rows_used=0, true_noise_used=False))
        return value
    return size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--stage', choices=('hard', 'both'), required=True)
    parser.add_argument('--depth-strata', type=int, choices=(1, 5), default=1)
    parser.add_argument('--continuous-state-guard', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    hook, runner, veto = [load(name) for name in ('gene_thinning_hook', 'membership_stage_panel', 'gene_veto_weights')]
    rec = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        hook=Path(hook.__file__), runner=Path(runner.__file__), veto=Path(veto.__file__),
        manifest=args.task.parent/'recovery_20261107/recovery_panel_v1/manifest.json'),
        dict(stage=args.stage, depth_strata=args.depth_strata, continuous_state_guard=args.continuous_state_guard,
             dispersion='reference-only common gene moment'),
        dict(scope='Non-oracle research recovery candidate; no calibration or biological validation'), mode='reference_gene_noise')
    estimates, splits, weights, argv = [], [], [], sys.argv
    (args.out/'folds').mkdir()
    (args.out/'weights').mkdir()
    def save_records():
        write_json(estimates, args.out/'estimates.json')
        write_json(splits, args.out/'split_records.json')
        write_json(weights, args.out/'weight_records.json')
    try:
        sys.argv = [str(runner.__file__), '--out', str(args.out/'panel'), '--stage', args.stage,
            '--panel-root', str(args.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(args.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        if args.continuous_state_guard:
            if args.stage != 'both':
                raise ValueError('Continuous-state guard requires the frozen both-stage candidate')
            sys.argv.append('--continuous-state-guard')
        context = veto.observe_weights(weights, args.out/'weights', apply=True, weight_units='effective') if args.stage == 'both' else nullcontext()
        with context:
            with hook.gene_native(splits, args.out/'folds', gene_size=estimator(estimates, args.depth_strata)):
                runner.main()
        summary = json.loads((args.out/'panel/summary.json').read_text())
        assert len(summary['cases']) == 5 and len(estimates) == 5
        write_json(dict(recovered=sum(row['original_assertion_passed'] for row in summary['cases']),
            cases=summary['cases'], fitted_reference_noises=estimates, production_changed=False,
            true_noise_used=False, scientific_criteria_changed=False), args.out/'report.json')
        save_records()
        finish_run(args.out, rec, backend='CPU reference-only dispersion and full gene-native workflow')
    except BaseException as error:
        save_records()
        finish_run(args.out, rec, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
