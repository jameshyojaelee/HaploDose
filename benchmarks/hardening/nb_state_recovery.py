"""Frozen original recovery with a count-native continuous-state safeguard."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('nb_guard_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    state, noise = load('nb_state_guard'), load('reference_gene_noise')
    records = []
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        guard=Path(state.__file__), noise=Path(noise.__file__),
        stage=Path(__file__).with_name('membership_weighted_stages.py'),
        runner=Path(__file__).with_name('membership_stage_panel.py'),
        thinning=Path(__file__).with_name('gene_thinning_hook.py'),
        veto=Path(__file__).with_name('gene_veto_weights.py')),
        dict(unchanged_support_thresholds=True, BIC_extra_parameters=2),
        dict(scope='Recovery prerequisite only; biological specimens=0'), mode='nb_state_recovery')
    original_load, argv = noise.load, sys.argv
    def with_guard(name):
        module = original_load(name)
        if name != 'membership_stage_panel':
            return module
        runner_load = module.load
        def patched_load(label, path):
            loaded = runner_load(label, path)
            if label == 'weighted_stages':
                original = loaded.weighted_stages
                @contextmanager
                def guarded(stage, trace, **kwargs):
                    assert kwargs.get('candidate_filter') is None
                    kwargs['candidate_filter'] = state.make_guard(records)
                    with original(stage, trace, **kwargs):
                        yield
                loaded.weighted_stages = guarded
            return loaded
        module.load = patched_load
        return module
    noise.load = with_guard
    try:
        sys.argv = [str(noise.__file__), '--out', str(args.out/'candidate'), '--task',
            str(args.task.parent/'depth_stratified_gene_noise_20261109'), '--stage', 'both', '--depth-strata', '5']
        noise.main()
        report = json.loads((args.out/'candidate/report.json').read_text())
        write_json(dict(recovered=report['recovered'], cases=report['cases'], guard_calls=len(records),
            false_split_criteria_changed=False, production_clone_changes=False), args.out/'report.json')
        write_json(records, args.out/'guard_records.json')
        finish_run(args.out, run, backend='CPU complete gene-count recovery and NB state comparison')
    except BaseException as error:
        write_json(records, args.out/'guard_records.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        noise.load, sys.argv = original_load, argv


if __name__ == '__main__':
    main()
