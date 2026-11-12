"""Bounded gene-native versus moment-bin continuous-state performance probe."""
import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import nbinom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('gene_state_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def observations(task, seed):
    prior = task.parent/'conditioned_gene_score_20261109/full_candidate_v2'
    trace = json.loads((prior/'stage_trace.json').read_text())
    selected = [(i, entry, call) for i, entry in enumerate(trace)
        if entry.get('case_seed') == seed and entry.get('accepted')
        and entry.get('source', '').startswith('event_chrom12_')
        for call in entry.get('calls', []) if call['context'] == 'depth_matched' and call.get('exact_segments')]
    assert len(selected) == 1
    index, entry, call = selected[0]
    lock = json.loads(Path(entry['proposal_lock']).read_text())
    with np.load(call['arrays']) as z:a = {k:z[k] for k in z.files}
    matches = []
    for record in json.loads((prior/'candidate/split_records.json').read_text()):
        if record['case_seed'] != seed:continue
        with np.load(record['path']) as z:
            if len(z['global_rows']) <= int(a['cells'].max()):continue
            if np.array_equal(z['global_rows'][a['cells']], a['global_rows']):
                matches.append(record)
    assert len(matches) == 1
    with np.load(matches[0]['path']) as z:
        prop = sparse.csr_matrix((z['prop_data'], z['prop_indices'], z['prop_indptr']), shape=z['prop_shape'])
        counts = prop[a['cells']][:, a['genes']].toarray()
        share = float(z['share'])
    use = np.asarray(lock['parameters']['selected'], int)
    np.testing.assert_array_equal(use, a['selected'])
    np.testing.assert_array_equal(lock['genes'], a['genes'])
    np.testing.assert_array_equal(lock['rows'], a['global_rows'])
    scale = share/(1-share)
    expected = a['expected']*scale*np.exp(np.asarray(lock['parameters']['offsets'])[:, None])
    sizes = np.broadcast_to(a['size']*scale, counts.shape)
    return counts[use], expected[use], sizes[use], a['q'][use] > .5, a['local_bin'], dict(
        source=entry['source'], trace_index=index, confirmation_counts_used=False,
        proposal_lock=entry['proposal_lock'], fold=matches[0]['path'],
        selected_global_rows=a['global_rows'][use].tolist())


def convolution_check():
    # Different success probabilities: moment matching is not NB closure.
    mean = np.array([.5, 4.]);size = np.array([1., 1.]);total = mean.sum()
    moment_size = total**2/np.sum(mean**2/size)
    exact = np.convolve(nbinom.pmf(np.arange(15), size[0], size[0]/(size[0]+mean[0])),
        nbinom.pmf(np.arange(15), size[1], size[1]/(size[1]+mean[1])))[:15]
    approximation = nbinom.pmf(np.arange(15), moment_size, moment_size/(moment_size+total))
    discrepancy = float(np.max(abs(exact-approximation)))
    operations = 128*15**2
    bound = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    assert discrepancy > bound
    return dict(moment_nb_is_not_exact=True, maximum_probability_discrepancy=discrepancy,
        rounding_budget=bound, independent_scipy_convolution=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--seed', type=int, choices=(20261031, 20261107, 20261108, 20261109, 20261110), required=True)
    parser.add_argument('--require-corrected-null', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    prior = args.task.parent/'conditioned_gene_score_20261109/full_candidate_v2'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        guard=Path(__file__).with_name('nb_state_guard.py'),
        numerical_reference=args.task.parent/'nb_state_recovery_20261109/numerics_v1/report.json',
        trace=prior/'stage_trace.json', splits=prior/'candidate/split_records.json')
    if args.require_corrected_null:
        inputs['counterexample'] = args.task/'boundary_v2/report.json'
        inputs['correction_design'] = args.task/'CORRECTION_DESIGN.md'
    rec = start_run(args.out, inputs, dict(seed=args.seed, orders=[24,48,96,192], BIC_extra_parameters=2),
        dict(scope='Conditional proposal comparison only; no full recovery/error or biological claim'), mode='gene_state_contrast')
    try:
        numerical = convolution_check()
        O, E, r, g, bins, provenance = observations(args.task, args.seed)
        np.savez_compressed(args.out/'proposal_observations.npz', counts=O, expected=E, size=r, labels=g, bins=bins)
        write_json(provenance, args.out/'OBSERVATIONS_LOCKED.json')
        records = []
        for arm in ('moment_bin', 'native_gene'):
            if arm == 'moment_bin':
                yy = np.stack([O[:,bins == j].sum(axis=1) for j in np.unique(bins)], axis=1)
                ee = np.stack([E[:,bins == j].sum(axis=1) for j in np.unique(bins)], axis=1)
                rr = ee**2/np.stack([(E[:,bins == j]**2/r[:,bins == j]).sum(axis=1) for j in np.unique(bins)], axis=1)
            else:yy, ee, rr = O, E, r
            studies = []
            for order in (24, 48, 96, 192):
                result = load('nb_state_guard').fit_comparison(yy, ee, rr, g, order)
                if args.require_corrected_null:
                    counter = next(x for x in json.loads(inputs['counterexample'].read_text())['records']
                        if x['seed'] == args.seed and x['arm'] == arm)
                    grid_value = counter['best_fixed_mean']['objectives'][(24,48,96).index(order)] if order <= 96 \
                        else counter['best_fixed_mean']['objectives'][-1]
                    assert result['null_objective'] <= grid_value+counter['rounding_budget']
                studies.append(result)
                write_json(dict(arm=arm, quadrature=studies), args.out/f'{arm}.json')
                if len(studies) >= 2:
                    discrepancy = abs(studies[-1]['BIC_advantage']-studies[-2]['BIC_advantage'])
                    if abs(result['BIC_advantage']) > 2*discrepancy:break
            resolved = bool(abs(result['BIC_advantage']) > 2*discrepancy)
            records.append(dict(arm=arm, cells=len(yy), features=yy.shape[1], quadrature=studies,
                resolved=resolved, passed=bool(resolved and result['BIC_advantage'] > 0)))
            write_json(dict(seed=args.seed, numerical=numerical, comparisons=records,
                original_threshold_changed=False, full_workflow_run=False, biological_specimens=0), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU gene-native versus moment-bin shared continuous-state likelihood')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':main()
