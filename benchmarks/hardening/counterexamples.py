"""Same independent counterexamples runnable against historical and candidate wheels."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import optimize, special, stats, sparse

import haplodose
from haplodose.bulk import reference_model
from haplodose.io import CountSource
from haplodose._native import bin_sums
from haplodose.kernels.estep import programme_scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    result = dict(version=haplodose.__version__, package=haplodose.__file__,
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    reset = np.zeros(40, bool)
    reset[[0, 20]] = True
    _, _, W = reference_model(np.tile(np.arange(1, 41, dtype=float), (4, 1)), 10, reset)
    result['constant_reference_rank'] = dict(observed=len(W), independent_expected=0)
    rng = np.random.default_rng(82)
    J, F = 30, 2
    W = rng.normal(0, .3, (F, J))
    log_rate = rng.normal(-np.log(J), .5, (3, J))
    phi = rng.uniform(.05, .6, J)
    R, prior = np.array([.15, .55, .30]), np.array([.8, 1.1])
    y = rng.negative_binomial(1/phi, 1/(1+phi*1000*special.softmax(log_rate[1]+np.array([.4,-.7])@W)))
    def objective(u):
        mean = y.sum()*special.softmax(log_rate+u@W, axis=1)
        return float(R @ stats.nbinom.logpmf(y[None], 1/phi, 1/(1+phi*mean)).sum(axis=1)
                     - .5*np.sum(prior*u*u))
    independent = optimize.minimize(lambda u: -objective(u), np.zeros(F), method='BFGS')
    score = programme_scores(y[None], np.array([y.sum()]), R[None], log_rate, phi, W, prior)[0]
    result['programme_objective'] = dict(independent_score=independent.x.tolist(), score=score.tolist(),
                                        independent_minus_package=objective(independent.x)-objective(score),
                                        arithmetic_budget=float((J+16)*np.finfo(np.float32).eps*(1+abs(objective(independent.x)))))
    for name, data, genes in [('fractional', [1.0001, 2.], ['A','B']),
                              ('duplicate_features', [1., 2.], ['A','A'])]:
        path = args.out/f'{name}.h5'
        with h5py.File(path,'x') as h:
            h['genes'], h['cell_ids'] = np.array(genes,dtype='S'), np.array(['C'],dtype='S')
            h['X/data'], h['X/indices'], h['X/indptr'] = data, [0,1], [0,2]
        try:
            with CountSource(path) as source:
                value = source.read(0,1).toarray().tolist()
            result[name] = dict(rejected=False, value=value, independent_expected='reject')
        except ValueError as error:
            result[name] = dict(rejected=True, error=str(error), independent_expected='reject')
    try:
        value, _ = bin_sums(sparse.csr_matrix([[2**24,1.]]), [0,0], 1)
        result['aggregate_precision'] = dict(rejected=False, value=float(value[0,0]), independent_expected=2**24+1)
    except ValueError as error:
        result['aggregate_precision'] = dict(rejected=True, error=str(error), independent_expected=2**24+1)
    (args.out/'results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
