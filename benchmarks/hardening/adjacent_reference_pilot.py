"""Fresh complete-workflow test of the frozen adjacent-reference noise method."""
import argparse
from contextlib import nullcontext
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


PHASE = 7
ARMS = ('production', 'adjacent_gene_score')


def load(name):
    spec = importlib.util.spec_from_file_location('conditioned_pilot_'+name,
        Path(__file__).with_name(name+'.py'))
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
    args = parser.parse_args()
    assert 0 <= args.start < args.stop <= 100
    args.out.mkdir(exist_ok=False)
    fresh, noise, hook, method, veto = [load(name) for name in
        ('gene_veto_fresh', 'adjacent_reference_noise', 'gene_thinning_hook',
         'gene_shared_confirmation', 'gene_veto_weights')]
    inputs = dict(script=Path(__file__), design=args.task/'PILOT_DESIGN.md',
        candidate_report=args.task/'candidate_v2/report.json',
        control_report=args.task/'strict_control_v2/report.json',
        frozen_manifest=args.task/'source_v2/source_manifest.json')
    for name in ('gene_veto_fresh', 'adjacent_reference_noise', 'gene_thinning_hook',
            'gene_shared_confirmation', 'gene_veto_weights', 'membership_weighted_stages',
            'gene_event_proposal', 'shared_membership_stages', 'shared_membership_hmm',
            'depth_conditioned_genes', 'conditioned_gene_score', 'predictive_gene_effect'):
        inputs[name] = Path(__file__).with_name(name+'.py')
    frozen = json.loads(inputs['frozen_manifest'].read_text())
    for name, path in inputs.items():
        if name in ('script', 'design', 'candidate_report', 'control_report', 'frozen_manifest'):
            continue
        assert hashlib.sha256(path.read_bytes()).hexdigest() == frozen['benchmarks/hardening/'+path.name]
    assert json.loads(inputs['candidate_report'].read_text())['recovered'] == 5
    assert json.loads(inputs['control_report'].read_text())['recovered'] == 5
    params = clones.CloneParams(split_mode='thin')
    rec = start_run(args.out, inputs, dict(phase=PHASE, condition=args.condition,
        start=args.start, stop=args.stop, arms=ARMS, params=asdict(params),
        reference_noise='disjoint_adjacent_depth_pairs', candidate_frozen=True),
        dict(scope='Fresh conditional complete-workflow test; biological n=0'),
        mode='adjacent_reference_pilot')
    results = []
    try:
        for rep in range(args.start, args.stop):
            seed, X, genes, metadata, positions, arrays = fresh.generate(rep, PHASE, args.condition)
            case = args.out/f'seed_{seed}'
            case.mkdir()
            fresh.write_input(case/'input', X, genes, metadata, positions)
            np.savez_compressed(case/'generation.npz', **arrays)
            for arm in ARMS:
                out = case/arm
                out.mkdir()
                for name in ('folds', 'weights', 'arrays'):
                    (out/name).mkdir()
                estimates, split_trace, stage_trace, veto_trace = [], [], [], []
                with veto.observe_weights(veto_trace, out/'weights', apply=True,
                        weight_units='effective') if arm != 'production' else nullcontext():
                    with hook.gene_native(split_trace, out/'folds',
                            gene_size=noise.estimator(estimates, 5), case_seeds=(seed,)) \
                            if arm != 'production' else nullcontext() as registry:
                        with method.confirmation(registry, 'depth_gene_score', stage_trace, out/'arrays') \
                                if arm != 'production' else nullcontext():
                            pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv',
                                case/'input/positions.tsv', out/'fit', params=params, chunk=97,
                                counts_kind='raw', genome_build='hg38',
                                coordinate_system='one-based-closed', gene_namespace='symbol')
                write_json(dict(estimates=estimates, splits=split_trace, stages=stage_trace,
                    veto=veto_trace), out/'trace.json')
                row = fresh.score(out/'fit', metadata, args.condition)
                row.update(replicate=rep, phase=PHASE, seed=seed, arm=arm, case=str(case),
                    fit_manifest=file_record(out/'fit/run_manifest.json'))
                results.append(row)
                write_json(results, args.out/'results.json')
                print(json.dumps(row), flush=True)
        finish_run(args.out, rec, backend='CPU two paired complete workflows, frozen adjacent-reference candidate')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
