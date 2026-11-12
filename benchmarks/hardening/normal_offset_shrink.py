"""Normal-only cross-validated anchor shrinkage; query and truth processed later."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def model_data(Y, L, base, nuisance, phi):
    return dict(Y=Y.astype(float), L=L, base=base, nuisance=nuisance, phi=phi,
        R=np.ones((len(Y), 1)), doses=np.zeros((1, len(base))), offset=np.zeros(len(base)))


def normal_selection(N, L, base, nuisance, phi, helper, estimator):
    """Only normal observations enter folds, selection and full normal refitting."""
    folds = [np.arange(parity, len(N), 2) for parity in (0, 1)]
    predictions, optimization = [], []
    for rows in folds:
        nb, pool, diagnostic = estimator.estimate(N[rows], L[rows], base, nuisance[rows], phi, helper)
        predictions.append(dict(NB=nb, pooled=pool))
        optimization.append(diagnostic)
    full_nb, full_pool, full_diagnostic = estimator.estimate(N, L, base, nuisance, phi, helper)
    result, arrays = {}, {}
    for family, full in (('NB', full_nb), ('pooled', full_pool)):
        validation = [model_data(N[rows], L[rows], base, nuisance[rows], phi) for rows in folds]

        def score(t):
            value, budget = 0., 0.
            for held in (0, 1):
                gain, _, error = helper.gain_gradient(validation[held], t*predictions[1-held][family])
                value += gain
                budget += error
            return value, budget

        eps = np.finfo(float).eps
        coarse = minimize_scalar(lambda t: -score(t)[0], bounds=(0., 1.), method='bounded',
            options=dict(xatol=eps**.5, maxiter=500))
        fine = minimize_scalar(lambda t: -score(t)[0], bounds=(0., 1.), method='bounded',
            options=dict(xatol=eps**.75, maxiter=500))
        if not coarse.success or not fine.success:
            raise ArithmeticError('Unresolved scalar shrinkage optimizer')
        candidates = [(0., *score(0.)), (1., *score(1.)), (float(fine.x), *score(fine.x))]
        t, value, budget = max(candidates, key=lambda v: v[1])
        coarse_value, coarse_budget = score(coarse.x)
        if value < coarse_value-budget-coarse_budget:
            raise ArithmeticError('Fine scalar selection loses the coarse validation objective')
        result[family] = dict(t=t, heldout_gain=value, numerical_budget=budget,
            endpoint_zero_gain=candidates[0][1], endpoint_one_gain=candidates[1][1],
            coarse_t=float(coarse.x), fine_t=float(fine.x), coarse_gain=coarse_value,
            normal_prediction_refinement=float(np.max(np.abs((coarse.x-fine.x)*full))),
            full_fit=full_diagnostic, half_fits=optimization,
            scope='Full normal refit is descriptive; heldout scores remain those of half fits')
        for label, value in (('full', full), ('even', predictions[0][family]), ('odd', predictions[1][family])):
            arrays[family+'_'+label] = value
    digest = hashlib.sha256()
    for value in (N, L, base, nuisance, phi):
        value = np.ascontiguousarray(value)
        digest.update(str((value.shape, value.dtype.str)).encode())
        digest.update(value.tobytes())
    result['normal_input_sha256'] = digest.hexdigest()
    return result, arrays


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--observations', type=Path)
    p.add_argument('--truth', type=Path)
    p.add_argument('--information-scaled', action='store_true')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    panel = a.task/'panel_v1'
    observations, truth = a.observations or panel/'observations.npz', a.truth or panel/'targets.npz'
    root = Path(__file__).parent
    helper, estimator = load('objective', root/'offset_nb_objective.py'), load('estimator', root/'reference_offset_panel.py')
    info = load('information', root/'normal_offset_information.py') if a.information_scaled else None
    rec = start_run(a.out, dict(script=Path(__file__), helper=root/'offset_nb_objective.py',
        estimator=root/'reference_offset_panel.py', design=a.task/'DESIGN.md',
        observations=observations, truth=truth, metadata=panel/'metadata.json'),
        dict(information_scaled=True) if a.information_scaled else {},
        dict(scope='Normal-only selection; known exposure/nuisance; biological n=0'), mode='normal_offset_shrink')
    try:
        observed = np.load(observations)
        keys = json.loads((panel/'metadata.json').read_text())['keys']
        references, reference_arrays = {}, {}
        for key in keys:
            values = {name: observed[key+'__'+name] for name in ('N', 'Ln', 'base', 'nn', 'phi')}
            selected, arrays = normal_selection(values['N'], values['Ln'], values['base'], values['nn'],
                values['phi'], helper, estimator)
            if info is not None:
                selected['information'] = info.selection(values['N'], values['Ln'], values['base'],
                    values['nn'], values['phi'], helper, arrays)
            references[key] = selected
            reference_arrays.update({key+'__'+name: value for name, value in arrays.items()})
        # Seal normal artifacts before query counts or scoring truth are accessed.
        (a.out/'references.json').write_text(json.dumps(references, indent=2, sort_keys=True)+'\n')
        np.savez_compressed(a.out/'reference_arrays.npz', **reference_arrays)
        target, fitted, records = np.load(truth), {}, []
        for key in keys:
            condition, rep = key.split('__')
            qnb, qpool, _ = estimator.estimate(observed[key+'__Q'], observed[key+'__Lq'], observed[key+'__base'],
                observed[key+'__nq'], observed[key+'__phi'], helper)
            normal = {family: reference_arrays[key+'__'+family+'_full'] for family in ('NB', 'pooled')}
            profiles = dict(NB_zero=qnb, NB_anchor=estimator.centre(qnb-normal['NB']),
                NB_shrink=estimator.centre(qnb-references[key]['NB']['t']*normal['NB']),
                pooled_zero=qpool, pooled_anchor=estimator.centre(qpool-normal['pooled']),
                pooled_shrink=estimator.centre(qpool-references[key]['pooled']['t']*normal['pooled']))
            if info is not None:
                profiles.update(NB_information=estimator.centre(qnb-references[key]['information']['NB']['t']*normal['NB']),
                    pooled_information=estimator.centre(qpool-references[key]['information']['pooled']['t']*normal['pooled']))
            for method, profile in profiles.items():
                fitted[key+'__'+method] = profile
                dna = target[key+'__DNA_linked_RNA']
                rna = target[key+('__RNA_zero' if method.endswith('zero') else '__RNA_anchor')]
                records.append(dict(condition=condition, replicate=int(rep), method=method,
                    DNA_linked_RNA_MSE=float(np.mean((profile-dna)**2)),
                    own_relative_RNA_MSE=float(np.mean((profile-rna)**2)),
                    region_contrast=float(profile[12:22].mean()-np.r_[profile[:12], profile[22:]].mean())))
                if info is not None:
                    truth_contrast = dna[12:22].mean()-np.r_[dna[:12], dna[22:]].mean()
                    records[-1]['regional_squared_error'] = float((records[-1]['region_contrast']-truth_contrast)**2)
        np.savez_compressed(a.out/'fitted.npz', **fitted)
        (a.out/'records.json').write_text(json.dumps(records, indent=2)+'\n')
        summary = []
        for condition in ('null', 'gain', 'shared', 'mismatched', 'state'):
            for method in profiles:
                rows = [v for v in records if v['condition'] == condition and v['method'] == method]
                summary.append(dict(condition=condition, method=method, simulation_draws=len(rows),
                    mean_DNA_linked_RNA_MSE=float(np.mean([v['DNA_linked_RNA_MSE'] for v in rows])),
                    mean_own_relative_RNA_MSE=float(np.mean([v['own_relative_RNA_MSE'] for v in rows])),
                    mean_region_contrast=float(np.mean([v['region_contrast'] for v in rows])),
                    positive_region_contrast=sum(v['region_contrast'] > 0 for v in rows)))
                if info is not None:
                    summary[-1]['mean_regional_squared_error'] = float(np.mean([v['regional_squared_error'] for v in rows]))
        (a.out/'report.json').write_text(json.dumps(dict(summary=summary, biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False,
            scope='Conditional simulations; no clone, biological or uncertainty guarantee'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU sealed normal-only selection then independent query fitting/scoring')


if __name__ == '__main__':
    main()
