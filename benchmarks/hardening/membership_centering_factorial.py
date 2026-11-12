"""Separate oracle membership information from component-rate centering.

Postfit mechanism diagnosis only. Neither simulated labels nor oracle rates enter
the package, proposal fitting, confirmation decisions or calibration experiments.
"""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import nbinom
from haplodose.clones import CloneParams, _thin_diff_params
from haplodose.profile import DIFF_LEVELS, DIFF_PRIOR, credible_difference_runs
from haplodose.provenance import start_run, finish_run, file_record


def reference(counts, exposure, size, q, pi, effects, prior, change, left, right):
    """Independent full-PMF sum over every membership and state on a small fixture."""
    cells, length = counts.shape
    states = list(itertools.product(range(3), repeat=length))
    transition = np.array([[1-change if i == j else change*prior[j]/(1-prior[i])
                            for j in range(3)] for i in range(3)])
    prior_mass, terms = [], []
    null = float(nbinom.logpmf(counts, size, size/(size+exposure)).sum())
    for path in states:
        mass = transition[1, path[0]] if left else prior[path[0]]
        for a, b in zip(path[:-1], path[1:]):
            mass *= transition[a, b]
        if right:
            mass *= transition[path[-1], 1]
        prior_mass.append(mass)
        likelihood = []
        for membership in itertools.product((0, 1), repeat=cells):
            weight = math.prod(q[c] if z else 1-q[c] for c, z in enumerate(membership))
            if weight == 0:
                continue
            means = np.array([[exposure[c, j]*math.exp(z*effects[s]) /
                               (1-pi[j]+pi[j]*math.exp(effects[s]))
                               for j, s in enumerate(path)] for c, z in enumerate(membership)])
            likelihood.append(math.log(weight)+float(nbinom.logpmf(counts, size, size/(size+means)).sum())-null)
        terms.append(math.log(mass)+float(logsumexp(likelihood)))
    evidence = float(logsumexp(terms)-math.log(math.fsum(prior_mass)))
    weights = np.exp(np.asarray(terms)-logsumexp(terms))
    posterior = np.array([[math.fsum(w for path, w in zip(states, weights) if path[j] == s)
                           for s in range(3)] for j in range(length)])
    # Accumulation scale includes full PMFs, paths, cells and bins. This is a
    # numerical allowance, never a posterior/support decision threshold.
    allowance = 64*np.finfo(float).eps*len(states)*cells*length*max(1., abs(null), max(abs(x) for x in terms))
    return evidence, posterior, allowance


def checks(model):
    counts = np.array([[2., 9., 3.], [8., 4., 1.], [1., 2., 6.]])
    exposure = np.array([[3., 5., 4.], [5., 3., 2.], [2., 4., 5.]])
    size = np.full_like(counts, 4.); effects = np.array([-.5, 0., .5])
    prior = np.array([.1, .8, .1]); soft = np.array([.2, .65, .4]); truth = np.array([0., 1., 0.])
    records = []
    for q in (soft, truth):
        for center_q in (soft, truth):
            pi = exposure.T@center_q/exposure.sum(axis=0)
            for left, right in itertools.product((False, True), repeat=2):
                fitted = model.build(counts, exposure, size, q, effects, prior, .07,
                    left_neutral=left, right_neutral=right, mixture_proportion=pi)
                actual = model.enumerate_states(fitted)
                evidence, posterior, allowance = reference(counts, exposure, size, q, pi, effects, prior, .07, left, right)
                np.testing.assert_allclose(actual['posterior'], posterior, rtol=0, atol=allowance)
                if abs(actual['log_evidence']-evidence) > allowance:
                    raise AssertionError('Independent full-PMF evidence disagrees')
                perm = np.array([2, 0, 1])
                permuted = model.enumerate_states(model.build(counts[perm], exposure[perm], size[perm], q[perm],
                    effects, prior, .07, left_neutral=left, right_neutral=right, mixture_proportion=pi))
                np.testing.assert_allclose(actual['posterior'], permuted['posterior'], rtol=0, atol=allowance)
                records.append(dict(left=left, right=right, allowance=allowance,
                    max_posterior_discrepancy=float(np.max(np.abs(actual['posterior']-posterior))),
                    evidence_discrepancy=abs(actual['log_evidence']-evidence)))
    return records


def main():
    p = argparse.ArgumentParser(); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--previous-oracle', type=Path, required=True); p.add_argument('--task', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('shared_membership_hmm.py')
    spec = importlib.util.spec_from_file_location('factorial_exact_model', path)
    model = importlib.util.module_from_spec(spec); spec.loader.exec_module(model)
    old_manifest = json.loads((a.previous_oracle/'run_manifest.json').read_text())
    if old_manifest['status'] != 'complete':
        raise ValueError('Incomplete original oracle investigation')
    inputs = {k: Path(r['path']) for k, r in old_manifest['input']['files'].items()}
    inputs.update(factorial=Path(__file__), current_model=path,
                  previous_manifest=a.previous_oracle/'run_manifest.json', previous_results=a.previous_oracle/'summary.json')
    for record in old_manifest['input']['files'].values():
        if file_record(record['path'])['sha256'] != record['sha256']:
            raise ValueError('Frozen oracle inputs changed')
    rec = start_run(a.out, inputs, dict(seeds=[20261108, 20261109, 20261110]),
        dict(scope='Postfit factorial attribution, no deployable oracle or biological claim',
             factors=['membership probabilities in marginal likelihood', 'exposure-weighted component-rate centering']),
        mode='membership_centering_factorial')
    try:
        (a.out/'numerical_checks.json').write_text(json.dumps(checks(model), indent=2))
        previous = json.loads(inputs['previous_results'].read_text()); rows = []
        pp = _thin_diff_params(CloneParams(split_mode='thin').profile)
        for seed in (20261108, 20261109, 20261110):
            meta = pd.read_csv(inputs[f'{seed}:truth'], sep='\t')
            bins = pd.read_csv(inputs[f'{seed}:fit/bins.tsv'], sep='\t'); chrom = bins.chromosome.to_numpy()
            beta = float(pd.read_csv(inputs[f'{seed}:fit/specimens.tsv'], sep='\t').set_index('specimen').loc['S1', 'beta'])
            events = [x for x in json.loads(inputs[f'{seed}:trace.json'].read_text()) if x['stage'] == 'event_fit']
            if len(events) != 1:
                raise ValueError('Locked representative event set changed')
            for event in events:
                with np.load(inputs[f'{seed}:event{event["index"]}']) as z:
                    O, E, r, q, cols = [z[k] for k in ('Ytest', 'E', 'r', 'q', 'cols')]
                    selected = meta.iloc[z['rows'][z['cells']]]
                if not selected.specimen.eq('S1').all() or selected.cell_id.duplicated().any():
                    raise ValueError('Cell correspondence failed')
                truth = selected.truth_clone.eq('clone2').to_numpy(float)
                probabilities = dict(soft=q, oracle_truth=truth)
                proportions = {name: np.clip(E.T@value/E.sum(axis=0), 0., 1.) for name, value in probabilities.items()}
                reference_likelihoods = {}
                for membership, probability in probabilities.items():
                    for centering, pi in proportions.items():
                        fitted = model.build(O, E, r, probability, DIFF_LEVELS*beta*np.log(1.5), DIFF_PRIOR, pp.p_change,
                            left_neutral=bool(cols[0] > 0 and chrom[cols[0]-1] == chrom[cols[0]]),
                            right_neutral=bool(cols[-1]+1 < len(chrom) and chrom[cols[-1]+1] == chrom[cols[-1]]),
                            mixture_proportion=pi)
                        if centering in reference_likelihoods:
                            np.testing.assert_array_equal(fitted['ll'], reference_likelihoods[centering])
                        else:
                            reference_likelihoods[centering] = fitted['ll']
                        result = model.enumerate_states(fitted)
                        reset = np.r_[True, np.zeros(len(cols)-1, bool)]
                        cores = [[int(cols[x]), int(cols[y]), int(s)] for x, y, s in
                                 credible_difference_runs(result['path'], result['posterior'], reset, pp)]
                        if membership == centering:
                            old = next(x for x in previous if x['seed'] == seed and x['membership'] == membership)
                            np.testing.assert_array_equal(result['posterior'], old['posterior'])
                            if cores != old['cores'] or result['log_evidence'] != old['log_evidence']:
                                raise AssertionError('Original oracle result not exactly reproduced')
                        row = dict(seed=seed, membership=membership, centering=centering,
                            log_evidence=result['log_evidence'], cores=cores, posterior=result['posterior'].tolist(),
                            mixture_proportion=pi.tolist(), cell_ids=selected.cell_id.tolist(),
                            membership_accuracy=float(np.mean((probability > .5) == truth)),
                            conditional_confirmation_accuracy=float(np.mean((result['membership'] > .5) == truth)),
                            minimum_non_neutral_posterior=float(np.min(1-result['posterior'][:, 1])))
                        rows.append(row)
                        print(json.dumps({k: v for k, v in row.items() if k not in ('posterior', 'mixture_proportion', 'cell_ids')}), flush=True)
        if len(rows) != 12:
            raise AssertionError('Factorial cases missing')
        (a.out/'summary.json').write_text(json.dumps(rows, indent=2))
        baseline = json.loads((a.task/'baseline/source_manifest.json').read_text()); core = 0
        for relative, digest in baseline.items():
            if file_record(a.task/'baseline'/relative)['sha256'] != digest:
                raise AssertionError('Baseline changed')
            if relative.startswith(('python/', 'tests/', 'crates/')):
                if file_record(a.task.parents[1]/relative)['sha256'] != digest:
                    raise AssertionError('Production/test/native files changed')
                core += 1
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256'] != record['sha256']:
                raise AssertionError('Frozen inputs changed')
        (a.out/'preservation.json').write_text(json.dumps(dict(baseline_files=len(baseline), unchanged_core_files=core,
            verified_inputs=len(rec['input']['files']), original_diagonal_results_exact=6, factorial_cells=12), indent=2))
        finish_run(a.out, rec, backend='cpu')
    except Exception as error:
        finish_run(a.out, rec, error=error); raise


if __name__ == '__main__':
    main()
