"""Known-generator count membership reference, never used by inference."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import linear_sum_assignment, minimize_scalar
from scipy.special import expit, logsumexp, roots_hermitenorm
from scipy.stats import nbinom

from haplodose.provenance import start_run, finish_run, file_record


DELTA = .6 * np.log(1.5)
SIGMA = .25


def log_likelihood(t, counts, depth, label, shape):
    """Four NB sums, dropping only class-independent combinatorial constants.

    t is standardized state. Within each quarter all NB probabilities and
    dispersions agree, so summing counts is sufficient conditional on true depth.
    Whole-genome rate renormalization affects all four quarters.
    """
    t = np.atleast_1d(t)
    amplitude = SIGMA * t + label * DELTA
    log_weights = np.stack((np.full_like(t, DELTA), np.zeros_like(t),
                            amplitude, np.zeros_like(t)), axis=-1)
    log_rate = log_weights - logsumexp(log_weights, axis=-1, keepdims=True)
    log_mean = np.log(depth)[:, None, None] + log_rate[None]
    normalizer = np.logaddexp(np.log(shape), log_mean)
    return np.sum(counts[:, None, :] * (log_mean - normalizer)
                  - shape * normalizer, axis=-1)


def evidence(counts, depth, shape, order):
    nodes, weights = roots_hermitenorm(order)
    return np.column_stack([logsumexp(log_likelihood(nodes, counts, depth, label, shape)
                             + np.log(weights)[None], axis=1) - .5*np.log(2*np.pi)
                            for label in (0, 1)])


def independent_integral(counts, depth, shape, label):
    # Adaptive integration of the full NB densities, written independently
    # of log_likelihood and the Gauss-Hermite nodes.
    def kernel(t):
        weights = np.array([np.exp(DELTA), 1., np.exp(SIGMA*t + label*DELTA), 1.])
        means = depth * weights / weights.sum()
        return float(nbinom.logpmf(counts, shape, shape/(shape+means)).sum() - t*t/2)
    peak = minimize_scalar(lambda t: -kernel(t), bounds=(-12., 12.), method='bounded')
    if not peak.success:
        raise ArithmeticError('Adaptive reference mode search failed')
    centre = kernel(peak.x)
    # Finite integration avoids overflow in the independent probability-scale
    # implementation. A separate Gaussian-tail/likelihood bound is recorded.
    value = 0.
    error = 0.
    for lo, hi in ((-12., peak.x), (peak.x, 12.)):
        integral, err = quad(lambda t: np.exp(kernel(t)-centre), lo, hi,
                             epsabs=2.**-40, epsrel=2.**-35, limit=200)
        value += integral
        error += err
    from scipy.special import log_ndtr
    tail_log_bound = np.log(2.) + log_ndtr(-12.)
    log_value = np.log(value)+centre-.5*np.log(2*np.pi)
    relative_tail_bound = float(np.exp(tail_log_bound-log_value))
    if error >= value or relative_tail_bound >= 1:
        raise ArithmeticError('Independent quadrature cannot resolve likelihood')
    log_error = float(-np.log1p(-error/value)-np.log1p(-relative_tail_bound))
    return log_value, log_error, relative_tail_bound


def correct(truth, labels, eligible):
    table = np.array([[np.count_nonzero((truth==t)&(labels==g)) for g in eligible]
                      for t in (0, 1)]).reshape(2, len(eligible))
    return int(table[linear_sum_assignment(-table)].sum()) if table.size else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    source = root/'dev/joint_fresh_null_20261108'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
                  prior_audit=source/'audit_v1/run_manifest.json')
    assert json.loads(inputs['prior_audit'].read_text())['status'] == 'complete'
    records = {}
    for arm in ('joint_proposal', 'joint_full'):
        records[arm] = []
        for start in (0, 50, 100, 150):
            directory = source/'chunks_v1'/arm/f'eval_two_{start}'
            path = directory/'run_manifest.json'
            manifest = json.loads(path.read_text())
            assert manifest['status'] == 'complete'
            for key in ('records.jsonl',):
                record = manifest['output']['files'][key]
                assert file_record(record['path'])['sha256'] == record['sha256']
            inputs[f'{arm}:{start}:manifest'] = path
            inputs[f'{arm}:{start}:records'] = directory/'records.jsonl'
            records[arm].extend(json.loads(s) for s in inputs[f'{arm}:{start}:records'].read_text().splitlines())
    generator_record = manifest['input']['files']['global_rank_generator.py']
    assert file_record(generator_record['path'])['sha256'] == generator_record['sha256']
    inputs['generator'] = Path(generator_record['path'])
    spec = importlib.util.spec_from_file_location('frozen_generator', inputs['generator'])
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    run = start_run(args.out, inputs, dict(cases=200, cells=600, bins=240,
                    orders=[16,32,64,128,256], delta=DELTA, state_sd=SIGMA),
                    dict(scope='Conditional oracle attribution, not inference or biological validation',
                         unavailable_inputs='True per-cell depth, generating model, event region and component effects',
                         priors='Equal marginal class probability; balanced sorting additionally uses known class sizes'),
                    mode='count_membership_reference')
    try:
        rows = []
        quadratures = []
        posterior = np.empty((200, 600))
        convergence = np.empty((200, 4))
        for i in range(200):
            key = [20261107, 120, 2, i]
            assert all(records[arm][i]['seed_key'] == key for arm in records)
            y, base, phi, truth, latent = generator.generate(
                np.random.default_rng(np.random.SeedSequence(key)), 'two_clones', 600, 240, return_latent=True)
            reset = np.zeros(240, bool)
            reset[::20] = True
            for name, value in dict(counts=y, reference=base, phi=phi, truth=truth,
                                    reset=reset, state_effect=latent['state_effect'], mean=latent['mean']).items():
                assert all(generator.array_hash(value) == records[arm][i]['input_sha256'][name] for arm in records), f'Regenerated {name} hash differs in case {i}'
            assert np.all(phi == .15) and np.all(base == 1/240)
            counts = y.reshape(600, 4, 60).sum(axis=2, dtype=np.float64)
            shape = 60/.15
            if i == 0:
                sufficiency = []
                for cell in (0, 599):
                    differences = []
                    for t in (-2., 0., 2.):
                        for label in (0, 1):
                            weights = np.array([np.exp(DELTA), 1., np.exp(SIGMA*t+label*DELTA), 1.])
                            means = latent['depth'][cell]*weights/weights.sum()
                            full = nbinom.logpmf(y[cell], 1/.15,
                                                1/(1+.15*np.repeat(means/60,60)))
                            pooled = nbinom.logpmf(counts[cell], shape, shape/(shape+means))
                            differences.append(float(full.sum()-pooled.sum()))
                    # This invariant requires only class/state-independent
                    # combinatorial factors. Allow accumulation error scaled
                    # by the actual PMFs and the 240 summands.
                    operations = 64*240
                    unit = np.finfo(float).eps/2
                    budget = operations*unit/(1-operations*unit)*(1+abs(full).sum()+abs(pooled).sum())
                    discrepancy = float(np.ptp(differences))
                    assert discrepancy <= budget, 'Quarter sums lost membership likelihood information'
                    sufficiency.append(dict(cell=cell, max_constant_difference=discrepancy,
                                            accumulation_allowance=float(budget)))
                (args.out/'sufficiency_checks.json').write_text(json.dumps(sufficiency,indent=2))
            levels = [evidence(counts, latent['depth'], shape, order) for order in (16,32,64,128,256)]
            odds = [level[:,1]-level[:,0] for level in levels]
            posterior[i] = expit(odds[-1])
            convergence[i] = [np.max(abs(after-before)) for before, after in zip(odds, odds[1:])]
            if i in (0,99,199):
                # Select from counts-derived ambiguity only, before reading truth
                # or any fitted cell labels for the numerical reference.
                for cell in sorted({0,599,int(np.argmin(abs(odds[-1])))}):
                    independent = [independent_integral(counts[cell], latent['depth'][cell], shape, label) for label in (0,1)]
                    # Add back the class-independent NB normalizing constants.
                    from scipy.special import gammaln
                    constant = float(np.sum(gammaln(counts[cell]+shape)-gammaln(shape)-gammaln(counts[cell]+1)+shape*np.log(shape)))
                    difference = np.array([s[0] for s in independent])-(levels[-1][cell]+constant)
                    quadratures.append(dict(replicate=i,cell=cell,log_evidence_difference=difference.tolist(),
                        adaptive_log_error=[s[1] for s in independent],
                        prior_tail_relative_bound=[s[2] for s in independent],
                        odds_difference=float(difference[1]-difference[0])))
            predicted = (odds[-1]>0).astype(int)
            balanced = np.zeros(600,int)
            balanced[np.argsort(odds[-1],kind='stable')[-300:]] = 1
            amplitude = latent['state_effect']+DELTA*truth
            row = dict(replicate=i, known_depth_count_correct=int(np.sum(predicted==truth)),
                       known_depth_balanced_correct=int(np.sum(balanced==truth)),
                       noiseless_amplitude_correct=int(np.sum((amplitude>DELTA/2)==truth)),
                       count_decisions_changed_128_to_256=int(np.sum((odds[-2]>0)!=(odds[-1]>0))),
                       max_log_odds_change_128_to_256=float(convergence[i,-1]))
            for arm in records:
                record = records[arm][i]['assignments']
                assert file_record(record['path'])['sha256'] == record['sha256']
                # Assignment hashes are recorded in the already inventoried
                # source journals; include them explicitly in the result as well.
                with np.load(record['path']) as fitted:
                    assert np.array_equal(truth, fitted['simulated_CNA_labels'])
                    labels = fitted['raw_MAP']
                    eligible = fitted['eligible_tumour_groups']
                    row[arm+'_correct'] = correct(truth, labels, eligible)
                    row[arm+'_oracle_coarsening_correct'] = int(sum(max(np.count_nonzero((truth==t)&(labels==g)) for t in (0,1)) for g in eligible))
                    row[arm+'_assignment_sha256'] = record['sha256']
            rows.append(row)
            if (i+1)%25 == 0:
                print(json.dumps(dict(completed=i+1,total=200)),flush=True)
        frame = pd.DataFrame(rows)
        frame.to_csv(args.out/'cases.tsv',sep='\t',index=False)
        np.savez_compressed(args.out/'posterior.npz',known_depth_count_probability=posterior,
                            max_log_odds_changes=convergence)
        (args.out/'quadrature.json').write_text(json.dumps(quadratures,indent=2))
        report = dict(cases=200,biological_specimens=0,oracle_only=True,
                      mean_accuracy={c:float(frame[c].mean()/600) for c in frame if c.endswith('_correct')},
                      decisions_changed_128_to_256=int(frame.count_decisions_changed_128_to_256.sum()),
                      max_log_odds_changes_by_order=convergence.max(axis=0).tolist(),
                      max_independent_log_odds_discrepancy=max(abs(s['odds_difference']) for s in quadratures),
                      numerical_scope='Quadrature convergence and selected independent integrations; not a formal global error enclosure',
                      decision='No inference/default/threshold/claim promoted; original errors and power remain unchanged')
        (args.out/'report.json').write_text(json.dumps(report,indent=2))
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
