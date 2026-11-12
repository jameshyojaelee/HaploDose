"""Original recovery under one fixed adjacent-depth reference noise estimator."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('paired_recovery_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--mode', choices=('control', 'weighted', 'count_gene'), default='control')
    parser.add_argument('--checks-only', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    method, runner = load('paired_reference_noise'), load('gene_shared_recovery')
    record = start_run(args.out, dict(script=Path(__file__), method=Path(method.__file__), runner=Path(runner.__file__),
        design=args.task/'DESIGN.md', noise=Path(__file__).with_name('reference_gene_noise.py'),
        shared=Path(__file__).with_name('gene_shared_confirmation.py')),
        dict(mode=args.mode), dict(scope='Conditional reference-only gene-noise performance intervention; biological n=0'),
        mode='paired_reference_noise')
    original_load, argv = runner.load, sys.argv
    try:
        checks = method.checks()
        if args.checks_only:
            write_json(checks, args.out/'report.json')
        else:
            def altered(name):
                module = original_load(name)
                if name == 'reference_gene_noise' and args.mode != 'control':
                    module.estimator = method.estimator
                return module
            runner.load = altered
            if args.mode == 'weighted':
                noise = load('reference_gene_noise')
                noise.estimator = method.estimator
                sys.argv = [str(noise.__file__), '--out', str(args.out/'run'), '--task',
                    str(args.task.parent/'depth_stratified_gene_noise_20261109'), '--stage', 'both', '--depth-strata', '5']
                noise.main()
            else:
                mode = 'count_gene' if args.mode == 'count_gene' else 'control'
                sys.argv = [str(runner.__file__), '--out', str(args.out/'run'), '--task', str(args.task), '--mode', mode]
                runner.main()
            report = json.loads((args.out/'run/report.json').read_text())
            write_json(dict(checks=checks, mode=args.mode, recovered=report['recovered'], cases=report['cases'],
                noise_source='fit_reference_only', larger_null_power_prerequisite='5/5 unchanged original recovery',
                production_changed=False, scientific_criteria_changed=False), args.out/'report.json')
        finish_run(args.out, record, backend='CPU adjacent-depth reference moments and complete native recovery')
    except BaseException as error:
        finish_run(args.out, record, error=error)
        raise
    finally:
        runner.load, sys.argv = original_load, argv


if __name__ == '__main__':
    main()
