"""Fixed original recovery with exact shared-membership confirmation."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('gene_shared_runner_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--mode', choices=('control', 'moment_bin', 'gene', 'count_bin', 'count_gene', 'depth_conditioned', 'depth_gene_score'), required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    (args.out/'arrays').mkdir()
    method, noise = load('gene_shared_confirmation'), load('reference_gene_noise')
    run = start_run(args.out, dict(script=Path(__file__), method=Path(method.__file__),
        noise=Path(noise.__file__), design=args.task/'DESIGN.md',
        engine=Path(__file__).with_name('shared_membership_hmm.py'),
        stage=Path(__file__).with_name('membership_weighted_stages.py'),
        gene=Path(__file__).with_name('gene_event_proposal.py'),
        matching=Path(__file__).with_name('shared_membership_stages.py'),
        runner=Path(__file__).with_name('membership_stage_panel.py'),
        thinning=Path(__file__).with_name('gene_thinning_hook.py'),
        veto=Path(__file__).with_name('gene_veto_weights.py'),
        conditioned_model=Path(__file__).with_name('depth_conditioned_genes.py'),
        conditioned_gene_score=Path(__file__).with_name('conditioned_gene_score.py')),
        dict(mode=args.mode, depth_strata=5, maximum_paths=1000000),
        dict(scope='Conditional research recovery; no independent confirmation or genomic calibration'),
        mode='gene_shared_recovery')
    original_load, argv = noise.load, sys.argv
    registry = None
    trace = []
    def with_shared(name):
        module = original_load(name)
        if name == 'gene_thinning_hook':
            original = module.gene_native
            @contextmanager
            def thinning(*a, **k):
                nonlocal registry
                with original(*a, **k) as selected:
                    registry = selected
                    yield selected
            module.gene_native = thinning
        if name == 'membership_stage_panel':
            original = module.load
            def altered(name, path):
                stage = original(name, path)
                if path.name == 'membership_weighted_stages.py':
                    @contextmanager
                    def shared(stage_name, records, **kwargs):
                        if stage_name != 'both' or kwargs.get('candidate_filter') is not None:
                            raise ValueError('Only the retained both-stage candidate is supported')
                        if registry is None:
                            raise RuntimeError('Missing gene-thinning registry')
                        with method.confirmation(registry, args.mode, trace, args.out/'arrays'):
                            yield
                    stage.weighted_stages = shared
                return stage
            module.load = altered
        return module
    noise.load = with_shared
    try:
        sys.argv = [str(noise.__file__), '--out', str(args.out/'candidate'), '--task',
            str(args.task.parent/'depth_stratified_gene_noise_20261109'), '--stage', 'both', '--depth-strata', '5']
        noise.main()
        report = json.loads((args.out/'candidate/report.json').read_text())
        if args.mode == 'control':
            prior = args.task.parent/'depth_stratified_gene_noise_20261109/both_v1'
            assert report['recovered'] == 4
            old_report = json.loads((prior/'report.json').read_text())
            assert [r['criteria'] for r in old_report['cases']] == [r['criteria'] for r in report['cases']]
            for row in report['cases']:
                before = prior/'panel'/f"seed_{row['seed']}"/'fit'
                after = args.out/'candidate/panel'/f"seed_{row['seed']}"/'fit'
                for path in before.iterdir():
                    if path.suffix == '.npy':
                        np.testing.assert_array_equal(np.load(path), np.load(after/path.name))
                    elif path.name.endswith(('.tsv', '.tsv.gz')):
                        pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'),
                            pd.read_csv(after/path.name, sep='\t'), check_exact=True)
            for path in (prior/'folds').glob('*.npz'):
                with np.load(path) as before, np.load(args.out/'candidate/folds'/path.name) as after:
                    assert before.files == after.files
                    for name in before.files:
                        np.testing.assert_array_equal(before[name], after[name])
        write_json(dict(mode=args.mode, recovered=report['recovered'], cases=report['cases'],
            conditional_comparisons=sum(len(t.get('calls', [])) for t in trace),
            production_changed=False, scientific_criteria_changed=False), args.out/'report.json')
        write_json(trace, args.out/'stage_trace.json')
        finish_run(args.out, run, backend='CPU exact conditional cell-shared membership')
    except BaseException as error:
        write_json(trace, args.out/'stage_trace.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        noise.load, sys.argv = original_load, argv


if __name__ == '__main__':
    main()
