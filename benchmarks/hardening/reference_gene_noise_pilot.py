"""Fresh paired full-workflow error/power pilot for estimated gene splitting."""
import argparse
from contextlib import nullcontext
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    import importlib.util
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('reference_pilot_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--condition', choices=('matched', 'state', 'two_clone'), required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--stop', type=int, required=True)
    parser.add_argument('--phase', type=int, choices=(2, 3, 4, 5), default=2)
    parser.add_argument('--depth-strata', type=int, choices=(1, 5), default=1)
    parser.add_argument('--expected-recovery', type=int, choices=(3, 4), default=3)
    parser.add_argument('--candidate-task', type=Path)
    parser.add_argument('--design-file', type=Path)
    parser.add_argument('--record-statistic', action='store_true')
    args = parser.parse_args()
    assert 0 <= args.start < args.stop <= 400
    args.out.mkdir(exist_ok=False)
    fresh, noise, hook, stage, veto = [load(name) for name in
        ('gene_veto_fresh', 'reference_gene_noise', 'gene_thinning_hook', 'membership_weighted_stages', 'gene_veto_weights')]
    params = clones.CloneParams(split_mode='thin')
    candidate_task = args.candidate_task or args.task
    inputs = dict(script=Path(__file__), design=args.design_file or args.task/'PILOT_DESIGN.md',
        candidate_report=candidate_task/'both_v1/report.json')
    gate = load('specimen_workflow_gate') if args.record_statistic else None
    if gate is not None:
        inputs['statistic_source'] = Path(gate.__file__)
    for module in (fresh, noise, hook, stage, veto):
        inputs[Path(module.__file__).stem] = Path(module.__file__)
    rec = start_run(args.out, inputs, dict(phase=args.phase, depth_strata=args.depth_strata,
        expected_recovery=args.expected_recovery, start=args.start, stop=args.stop,
        condition=args.condition, params=asdict(params)),
        dict(scope='Fresh conditional full-workflow pilot; no biological guarantee'), mode='reference_gene_noise_pilot')
    results = []
    try:
        assert json.loads(inputs['candidate_report'].read_text())['recovered'] == args.expected_recovery
        for rep in range(args.start, args.stop):
            seed, X, genes, metadata, positions, arrays = fresh.generate(rep, args.phase, args.condition)
            case = args.out/f'seed_{seed}'
            case.mkdir()
            fresh.write_input(case/'input', X, genes, metadata, positions)
            np.savez_compressed(case/'generation.npz', **arrays)
            for arm in ('production', 'candidate'):
                out = case/arm
                out.mkdir()
                folds, weights = out/'folds', out/'weights'
                folds.mkdir()
                weights.mkdir()
                estimates, split_trace, stage_trace, veto_trace = [], [], [], []
                with hook.gene_native(split_trace, folds, gene_size=noise.estimator(estimates, args.depth_strata), case_seeds=(seed,)) if arm == 'candidate' else nullcontext():
                    with stage.weighted_stages('both', stage_trace) if arm == 'candidate' else nullcontext():
                        with veto.observe_weights(veto_trace, weights, apply=True, weight_units='effective') if arm == 'candidate' else nullcontext():
                            pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv', case/'input/positions.tsv',
                                out/'fit', params=params, chunk=97, counts_kind='raw', genome_build='hg38',
                                coordinate_system='one-based-closed', gene_namespace='symbol')
                write_json(dict(estimates=estimates, splits=split_trace, stages=stage_trace, veto=veto_trace), out/'trace.json')
                row = fresh.score(out/'fit', metadata, args.condition)
                if gate is not None:
                    row['raw_max_pair_score'] = gate.statistic(out/'fit')
                row.update(replicate=rep, seed=seed, arm=arm, case=str(case))
                results.append(row)
                write_json(results, args.out/'results.json')
                print(json.dumps(row), flush=True)
        finish_run(args.out, rec, backend='CPU paired complete production/candidate gene-count workflows')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
