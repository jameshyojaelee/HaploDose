"""Audit original recovery plus gene/bin/veto split identities and native RNG draws."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from haplodose.clones import thin_counts
from haplodose.io import CountSource
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def csr(saved, prefix):
    return sparse.csr_matrix((saved[prefix+'_data'], saved[prefix+'_indices'], saved[prefix+'_indptr']),
        shape=tuple(saved[prefix+'_shape']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    recovery = load('original_criteria', root/'offset_objective_audit.py')
    inputs = dict(script=Path(__file__), scorer=root/'offset_objective_audit.py', design=a.task/'DESIGN.md',
        gate=a.task/'law_v1/report.json', control=a.task/'control_v1/panel/summary.json',
        gene=a.task/'gene_v1/panel/summary.json', gene_trace=a.task/'gene_v1/trace.json',
        producer_manifest=a.task/'gene_v1/run_manifest.json',
        old_control=a.task.parent/'profile_transition_20261109/legacy_v1/panel/summary.json')
    for row in json.loads(inputs['gene_trace'].read_text()):
        inputs[f"fold:{row['model']}"] = Path(row['path'])
    for seed in (20261031, 20261107, 20261108, 20261109, 20261110):
        inputs[f'{seed}:counts'] = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/counts.h5'
        inputs[f'{seed}:truth'] = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv'
        inputs[f'{seed}:mapping'] = a.task/f'gene_v1/panel/seed_{seed}/fit/genes.tsv.gz'
        for arm in ('control', 'gene'):
            for name in ('cells.tsv.gz', 'differences.tsv', 'qc.json'):
                inputs[f'{seed}:{arm}:{name}'] = a.task/f'{arm}_v1/panel/seed_{seed}/fit'/name
    rec = start_run(a.out, inputs, {}, dict(scope='Known-noise mechanism and frozen original criteria, biological n=0'),
        mode='gene_thinning_audit')
    try:
        law = json.loads(inputs['gate'].read_text())
        assert law['all_checks_passed']
        example = law['exact_counterexample']
        joint = [[Fraction(v) for v in row] for row in example['low_joint']]
        assert joint[0][0]*joint[1][1]-joint[0][1]*joint[1][0] == Fraction(example['determinant']) != 0
        assert Fraction(example['covariance']) == 0
        control, control_outcomes = recovery.recovery(a.task, 'control', 'v1')
        gene, gene_outcomes = recovery.recovery(a.task, 'gene', 'v1')
        old = json.loads(inputs['old_control'].read_text())
        assert [(v['seed'], v['criteria'], v['original_assertion_passed']) for v in control['cases']] == [
            (v['seed'], v['criteria'], v['original_assertion_passed']) for v in old['cases']]
        traces = json.loads(inputs['gene_trace'].read_text())
        assert len(traces) == 15
        audited, source_cache = [], {}
        for row in traces:
            saved = np.load(row['path'])
            X, prop, veto = (csr(saved, name) for name in ('source', 'prop', 'veto'))
            seed = row['case_seed']
            if seed not in source_cache:
                with CountSource(inputs[f'{seed}:counts']) as source:
                    source_cache[seed] = sparse.vstack([chunk for _, _, chunk in source.chunks(97)]).tocsr()
            mapping = pd.read_csv(inputs[f'{seed}:mapping'], sep='\t')
            keys, gbin = saved['global_rows'], saved['gbin']
            assert np.array_equal(mapping.bin, gbin)
            original = source_cache[seed][keys][:, mapping.column.to_numpy(int)]
            assert (original-X).nnz == 0
            full, half = X.toarray(), prop.toarray()
            assert np.all(half >= 0) and np.all(half <= full)
            size, share, split_seed = float(saved['gene_size']), float(saved['share']), int(saved['seed'])
            assert size == row['oracle_gene_size'] == 10.
            assert share == row['share'] and split_seed == row['seed']
            # Independent replay of stream identity and ordering, not an
            # independence proof; that proof is the separate joint-law gate.
            replay = np.zeros_like(half)
            base = (split_seed % 2**64) << 64
            for i, key in enumerate(keys):
                nonzero = np.flatnonzero(full[i] > 0)
                rng = np.random.Generator(np.random.Philox(key=base+int(key)))
                probability = rng.beta(share*size, (1-share)*size, len(nonzero))
                replay[i, nonzero] = rng.binomial(full[i, nonzero].astype(np.int64), probability)
            assert np.array_equal(replay, half)
            aggregation = sparse.csr_matrix((np.ones(len(gbin)), (np.arange(len(gbin)), gbin)),
                shape=(len(gbin), saved['full_bins'].shape[1]))
            assert np.array_equal((X@aggregation).toarray(), saved['full_bins'])
            assert np.array_equal((prop@aggregation).toarray(), saved['prop_bins'])
            local = saved['veto_local']
            assert len(np.unique(local)) == len(local)
            assert (veto-(X-prop)[local]).nnz == 0
            assert np.array_equal((veto@aggregation).toarray(),
                saved['full_bins'][local]-saved['prop_bins'][local])
            selected = np.unique([0, len(keys)//2, len(keys)-1])
            permuted = selected[::-1]
            actual = thin_counts(full[permuted], keys[permuted], np.full(X.shape[1], size), split_seed, share)
            assert np.array_equal(actual, half[permuted])
            chunked = np.concatenate([thin_counts(full[part], keys[part], np.full(X.shape[1], size), split_seed, share)
                for part in np.array_split(selected, 2) if len(part)])
            assert np.array_equal(chunked, half[selected])
            assert row['fold_calls'] > 0 and row['veto_calls'] == 1
            audited.append(dict(model=row['model'], case_seed=seed, cells=X.shape[0], genes=X.shape[1],
                all_gene_draws_replayed=True, gene_bin_veto_conservation=True,
                subset_order_chunk_invariance=True, fold_calls=row['fold_calls']))
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, original_criteria_checked=60,
            control_replayed_exactly=True, recovery=control_outcomes+gene_outcomes, folds=audited,
            biological_specimens=0, production_changed=False, scientific_criteria_changed=False,
            scope='Known true gene noise; full-workflow recovery diagnostic only, no biological/error/power claim'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent raw/RNG replay, conservation and unchanged60 criteria')


if __name__ == '__main__':
    main()
