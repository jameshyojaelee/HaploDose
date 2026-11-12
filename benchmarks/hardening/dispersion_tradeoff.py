"""Fresh paired complete-binned-workflow test of latent-weight dispersion."""
import argparse
from contextlib import nullcontext
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import jax
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
from haplodose.clones import CloneParams, Rows, fit_specimen
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--arm', choices=['MAP', 'soft'], required=True)
    p.add_argument('--condition', choices=['matched_nb', 'correlated_state', 'two_clones'], required=True)
    p.add_argument('--start', type=int, choices=[0, 100], required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    if jax.default_backend() != 'cpu':
        raise RuntimeError('CPU backend required inside inference')
    rank_source = a.task.parent/'recovery_next_20261107/rank_source_v1/benchmarks/hardening'
    generator_path = rank_source/'global_rank_generator.py'
    assert file_record(generator_path)['sha256'] == '50f189829dedbe3d4018045f3dd954eeb0142b9485b6d24b796824a244d78fc0'
    generator = load('generator', generator_path)
    rank = load('rank', rank_source/'global_rank_gate.py')
    hook = load('dispersion', Path(__file__).with_name('dispersion_weighted.py'))
    params = CloneParams(n_programmes=0, spatial_neighbors=0, seed=20261107, split_mode='thin')
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'TRADEOFF_DESIGN.md',
        generator=generator_path, rank=rank_source/'global_rank_gate.py', hook=Path(hook.__file__),
        recovery=a.task/'recovery_MAP_v1/summary.json', exact_audit=a.task/'matched_audit_v1/run_manifest.json'),
        dict(arm=a.arm, condition=a.condition, start=a.start, stop=a.start+100, params=asdict(params),
             seed_key_format=[20261109, 501, 'condition_index', 'replicate']),
        dict(scope='Fresh known-reference binned workflow, biological n=0; ten new bounds, total family163'),
        mode='dispersion_tradeoff')
    try:
        control = json.loads((a.task/'recovery_MAP_v1/summary.json').read_text())
        assert next(v for v in control['cases'] if v['seed'] == 20261031)['original_assertion_passed']
        assert json.loads((a.task/'matched_audit_v1/run_manifest.json').read_text())['status'] == 'complete'
        for i in range(a.start, a.start+100):
            start = time.monotonic()
            key = [20261109, 501, generator.CONDITIONS.index(a.condition), i]
            rng = np.random.default_rng(np.random.SeedSequence(key))
            Y, base, phi, truth, latent = generator.generate(rng, a.condition, 600, 240, return_latent=True)
            reset = np.zeros(240, bool)
            reset[::20] = True
            hashes = {name: generator.array_hash(value) for name, value in dict(counts=Y, reference=base,
                phi=phi, truth=truth, mean=latent['mean'], state=latent['state_effect'], reset=reset).items()}
            trace = []
            with (hook.weighted_dispersion(trace, save_moments=False) if a.arm == 'soft' else nullcontext()):
                result = fit_specimen(Rows(Y, np.arange(len(Y))), Y.sum(axis=1, dtype=np.float64),
                    base[:, None], reset, phi, params)
                score, labels, eligible, differences = rank.final_difference_score(result, params.min_cells)
            table = np.array([[np.sum((truth == k) & (labels == group)) for group in eligible] for k in (0, 1)])
            if eligible:
                row, col = linear_sum_assignment(-table)
                accuracy = float(table[row, col].sum()/len(Y))
            else:
                accuracy = 0.
            path = a.out/f'case_{i:03d}.npz'
            np.savez_compressed(path, truth=truth, raw_MAP=labels, raw_R=result['R'], eligible=np.array(eligible, int),
                final_phi=result['model'].phi)
            record = dict(arm=a.arm, condition=a.condition, replicate=i, seed_key=key, input_sha256=hashes,
                raw_candidate_count=len(eligible), multigroup=len(eligible) > 1, exactly_two=len(eligible) == 2,
                all_cell_one_to_one_accuracy=accuracy, correct_two=len(eligible) == 2 and accuracy > .65,
                ARI=float(adjusted_rand_score(truth, labels)), final_difference_score=score,
                final_pair_diagnostics=differences, assignments=file_record(path),
                dispersion_trace=trace, seconds=time.monotonic()-start)
            with (a.out/'records.jsonl').open('a') as handle:
                handle.write(json.dumps(record)+'\n')
            print(json.dumps(dict(arm=a.arm, condition=a.condition, completed=i-a.start+1, total=100)), flush=True)
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU verified inside complete binned inference')


if __name__ == '__main__':
    main()
