"""Research shared-state prediction of confirmation from discovery counts."""
import argparse
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.special import expit, logsumexp, roots_hermitenorm
from scipy.stats import nbinom, norm

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load_guard():
    spec = importlib.util.spec_from_file_location('predictive_nb_state', Path(__file__).with_name('nb_state_guard.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marginal(y, expected, size, parameters, mixture, order):
    """Per-cell integrated NB log likelihood ratio relative to zero effect."""
    y, expected, size = np.broadcast_arrays(np.asarray(y, float), np.asarray(expected, float),
                                           np.asarray(size, float))
    parameters = np.asarray(parameters, float)
    means = parameters[:2] if mixture else parameters[:1]
    sigma = parameters[-1]
    if sigma < 0 or not all(np.isfinite(a).all() for a in (y, expected, size, parameters)):
        raise ValueError('Finite arrays and nonnegative latent SD required')
    if y.ndim != 2 or np.any(y < 0) or np.any(expected <= 0) or np.any(size <= 0):
        raise ValueError('Aligned counts and positive means/sizes required')
    if sigma == 0:
        nodes, weights = np.zeros(1), np.zeros(1)
    else:
        nodes, mass = roots_hermitenorm(order)
        weights = np.log(mass)-.5*np.log(2*np.pi)
    loge, logr = np.log(expected), np.log(size)
    base = np.logaddexp(logr, loge)
    components = []
    for mean in means:
        eta = mean+sigma*nodes
        term = y[:, :, None]*eta-(y+size)[:, :, None]*(
            np.logaddexp(logr[:, :, None], loge[:, :, None]+eta)-base[:, :, None])
        components.append(logsumexp(term.sum(axis=1)+weights, axis=1))
    joint = np.array(components).T
    if mixture:
        joint += np.array([-np.logaddexp(0., parameters[2]), -np.logaddexp(0., -parameters[2])])
    return logsumexp(joint, axis=1), np.exp(joint-logsumexp(joint, axis=1)[:, None])


def predictive(D, C, ED, EC, rD, rC, null, alternative, order):
    """Shared state in both folds: joint density divided by discovery density."""
    y, e, r = np.c_[D, C], np.c_[ED, EC], np.c_[rD, rC]
    nd, _ = marginal(D, ED, rD, null, False, order)
    nj, _ = marginal(y, e, r, null, False, order)
    ad, qd = marginal(D, ED, rD, alternative, True, order)
    aj, qj = marginal(y, e, r, alternative, True, order)
    score = (aj-ad)-(nj-nd)
    operations = 128*y.size*max(1, order)
    eps = np.finfo(float).eps
    magnitude = max(1., float(np.max(abs(np.r_[nd, nj, ad, aj]))))
    rounding = operations*eps/(1-operations*eps)*magnitude*len(D)
    def entropy(q):
        return float(np.mean(-np.sum(q*np.log(np.maximum(q, np.finfo(float).tiny)), axis=1)))
    return dict(score=float(math.fsum(map(float, score))), rounding_budget=float(rounding),
        per_cell=score, mixture_membership_discovery=qd[:, 1], mixture_membership_joint=qj[:, 1],
        discovery_membership_entropy=entropy(qd), joint_membership_entropy=entropy(qj))


def numerical_checks():
    D, C = np.array([[1., 2.], [3., 0.]]), np.array([[2., 1.], [0., 2.]])
    ED, EC = np.array([[1., 2.], [2., 1.]]), np.array([[2., 1.], [1., 2.]])
    rD, rC = np.full(D.shape, 2.), np.full(C.shape, 3.)
    null, alt = [0., 0.], [0., np.log(2.), 0., 0.]
    result = predictive(D, C, ED, EC, rD, rC, null, alt, 48)
    # Integer NB probabilities, exact binary component multiplier2 and pi1/2.
    def pmf(y, e, r, multiplier):
        p = Fraction(r, 1)/(Fraction(r, 1)+Fraction(e, 1)*multiplier)
        return Fraction(math.comb(y+r-1, y), 1)*p**r*(1-p)**y
    exact = []
    for i in range(len(D)):
        base_d = base_j = Fraction(1)
        mixed_d = mixed_j = Fraction(0)
        for j in range(D.shape[1]):
            base_d *= pmf(int(D[i,j]), int(ED[i,j]), 2, Fraction(1))
            base_j *= pmf(int(D[i,j]), int(ED[i,j]), 2, Fraction(1))*pmf(int(C[i,j]), int(EC[i,j]), 3, Fraction(1))
        for multiplier in (Fraction(1), Fraction(2)):
            md = mj = Fraction(1)
            for j in range(D.shape[1]):
                pd = pmf(int(D[i,j]), int(ED[i,j]), 2, multiplier)
                md *= pd
                mj *= pd*pmf(int(C[i,j]), int(EC[i,j]), 3, multiplier)
            mixed_d += md/2
            mixed_j += mj/2
        exact.append(math.log(float(mixed_j/mixed_d/(base_j/base_d))))
    exact_error = abs(result['score']-math.fsum(exact))
    assert exact_error <= result['rounding_budget']
    # Independent adaptive integration of actual PMFs, including their constants.
    null, alt = [.1, .17], [-.2, .3, -.4, .12]
    actual = predictive(D, C, ED, EC, rD, rC, null, alt, 96)
    independent, error_budget = [], 0.
    for i in range(len(D)):
        densities = []
        for parameters, mixture in ((null, False), (alt, True)):
            means = parameters[:2] if mixture else parameters[:1]
            probs = [1-expit(parameters[2]), expit(parameters[2])] if mixture else [1.]
            for joint in (False, True):
                value = err = 0.
                for mean, probability in zip(means, probs):
                    def integrand(z):
                        multiplier = np.exp(mean+parameters[-1]*z)
                        logp = nbinom.logpmf(D[i], rD[i], rD[i]/(rD[i]+ED[i]*multiplier)).sum()
                        if joint:
                            logp += nbinom.logpmf(C[i], rC[i], rC[i]/(rC[i]+EC[i]*multiplier)).sum()
                        return float(np.exp(logp)*norm.pdf(z))
                    part, er = quad(integrand, -12., 12., epsabs=256*np.finfo(float).eps,
                                    epsrel=256*np.finfo(float).eps, limit=200)
                    value += probability*part
                    err += probability*(er+2*norm.sf(12.))
                densities.append(value)
                error_budget += err/(value-err)
        independent.append(np.log(densities[3]/densities[2])-np.log(densities[1]/densities[0]))
    integration_error = abs(actual['score']-math.fsum(independent))
    budget = error_budget+actual['rounding_budget']
    assert integration_error <= budget
    swapped = [alt[1], alt[0], -alt[2], alt[3]]
    permutation = predictive(D, C, ED, EC, rD, rC, null, swapped, 96)
    assert abs(permutation['score']-actual['score']) <= actual['rounding_budget']+permutation['rounding_budget']
    return dict(exact_fraction_score=math.fsum(exact), exact_score_error=exact_error,
        exact_rounding_budget=result['rounding_budget'], independent_quad_error=integration_error,
        independent_quad_and_rounding_budget=budget, component_swap=True,
        integration_tail='12 standard deviations; PMF<=1 and explicit2*normal survival bound')


def read_proposal(case):
    trace = json.loads(Path(case['trace']).read_text())
    if isinstance(trace, dict):
        splits, trace = trace['splits'], trace['stages']
    else:
        splits = json.loads(Path(case['splits']).read_text())
    selected = [(i, t, c) for i, t in enumerate(trace)
        if t.get('case_seed') == case['seed'] and t.get('accepted')
        and t.get('source', '').startswith('event_chrom12_')
        for c in t.get('calls', []) if c['context'] == 'depth_matched' and c.get('exact_segments')]
    assert len(selected) == 1
    index, entry, call = selected[0]
    lock = json.loads(Path(entry['proposal_lock']).read_text())
    fields = ('cells', 'global_rows', 'genes', 'expected', 'size', 'selected', 'q', 'local_bin')
    with np.load(call['arrays']) as z:
        a = {key:z[key] for key in fields}  # Confirmation counts not opened.
    matches = []
    for record in splits:
        if record['case_seed'] != case['seed']:
            continue
        with np.load(record['path']) as z:
            if (len(a['cells']) and a['cells'].min() >= 0 and a['cells'].max() < len(z['global_rows'])
                    and np.array_equal(z['global_rows'][a['cells']], a['global_rows'])):
                matches.append(record)
    assert len(matches) == 1
    with np.load(matches[0]['path']) as z:
        proposal = sparse.csr_matrix((z['prop_data'], z['prop_indices'], z['prop_indptr']), shape=z['prop_shape'])
        D = proposal[a['cells']][:, a['genes']].toarray()
        share = float(z['share'])
    np.testing.assert_array_equal(a['selected'], lock['parameters']['selected'])
    np.testing.assert_array_equal(a['global_rows'], lock['rows'])
    np.testing.assert_array_equal(a['genes'], lock['genes'])
    use = a['selected']
    EC = a['expected']*np.exp(np.asarray(lock['parameters']['offsets'])[:, None])
    rC = np.broadcast_to(a['size'], EC.shape)
    return D[use], EC[use]*share/(1-share), rC[use]*share/(1-share), a['q'][use] > .5, a, dict(
        case=case, source=entry['source'], trace_index=index, fold=matches[0]['path'],
        confirmation=call['arrays'], proposal_lock=entry['proposal_lock'], selected=use.tolist(),
        selected_global_rows=a['global_rows'][use].tolist()), EC[use], rC[use]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--case', type=int)
    parser.add_argument('--numerics', action='store_true')
    parser.add_argument('--summary', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', cases=args.task/'CASES.json',
                  fitting=Path(__file__).with_name('nb_state_guard.py'))
    prepared = None
    if args.case is not None:
        case = json.loads(inputs['cases'].read_text())[args.case]
        prepared = read_proposal(case)
        provenance = prepared[5]
        inputs.update({key:Path(provenance[key]) for key in ('fold', 'confirmation', 'proposal_lock')})
        inputs['trace'] = Path(case['trace'])
        if 'splits' in case:
            inputs['splits'] = Path(case['splits'])
        inputs['numerical_checks'] = args.task/'numerics_v1/report.json'
    if args.summary:
        for i in range(6):
            inputs[f'case{i}'] = args.task/f'case_{i:02d}_v1/report.json'
    rec = start_run(args.out, inputs, dict(case=args.case, orders=[24,48,96]),
        dict(scope='Research predictive safeguard screen; biological n=0; no5% test or full-workflow guarantee'),
        mode='state_predictive_confirmation')
    try:
        if args.numerics:
            write_json(numerical_checks(), args.out/'report.json')
        elif args.summary:
            cases = json.loads(inputs['cases'].read_text())
            results = [json.loads((args.task/f'case_{i:02d}_v1/report.json').read_text()) for i in range(len(cases))]
            gate = all(r['resolved'] and ((r['score'] > 0) if c['kind']=='original_recovery' else (r['score'] < 0))
                       for c, r in zip(cases, results))
            write_json(dict(cases=results, screen_passed=gate, full_workflow_run=False,
                independent_evaluation=False, false_split_control_established=False,
                original_tests_or_scientific_criteria_changed=False,
                decision='Continue original full recovery and fresh null/power only' if gate else
                         'Reject predictive safeguard expansion; do not tune threshold'), args.out/'report.json')
        else:
            cases = json.loads(inputs['cases'].read_text())
            case = cases[args.case]
            D, ED, rD, g, metadata, provenance, EC, rC = prepared
            guard = load_guard()
            models = [guard.fit_comparison(D, ED, rD, g, order) for order in (24,48,96)]
            write_json(dict(case=case, models=models, provenance=provenance,
                confirmation_counts_opened=False, proposal_identity=file_record(Path(provenance['fold']))),
                args.out/'DISCOVERY_MODELS_LOCKED.json')
            with np.load(provenance['confirmation']) as z:
                C = z['counts'][metadata['selected']]
            with np.load(provenance['fold']) as z:
                full = sparse.csr_matrix((z['source_data'], z['source_indices'], z['source_indptr']), shape=z['source_shape'])
            total = full[metadata['cells']][:, metadata['genes']].toarray()[metadata['selected']]
            np.testing.assert_array_equal(total, D+C)
            rows = []
            for model in models:
                result = predictive(D, C, ED, EC, rD, rC, model['null'], model['alternative'], model['order'])
                arrays = {key:result.pop(key) for key in ('per_cell','mixture_membership_discovery','mixture_membership_joint')}
                np.savez_compressed(args.out/f"predictive_{model['order']}.npz", **arrays)
                rows.append(dict(order=model['order'], discovery_BIC=model['BIC_advantage'], **result))
            discrepancy = max(abs(rows[j]['score']-rows[j-1]['score']) for j in (1,2))
            uncertainty = 2*discrepancy+sum(r['rounding_budget'] for r in rows)
            resolved = bool(abs(rows[-1]['score']) > uncertainty)
            write_json(dict(case=case, source=provenance['source'], cells=len(D), genes=D.shape[1],
                quadrature=rows, score=rows[-1]['score'], numerical_uncertainty_diagnostic=uncertainty,
                resolved=resolved, convergence_is_not_formal_bound=True,
                proposed_predictive_guard_retains=bool(resolved and rows[-1]['score'] > 0),
                fit_reads_confirmation=False, fold_count_conservation_exact=True,
                full_workflow_run=False, false_split_control_established=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU shared-state NB discovery fit and joint predictive integration')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
