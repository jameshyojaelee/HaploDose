"""Fixed known-gene-noise thinning intervention, unchanged original five criteria."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import sys

from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--mode', choices=['control', 'gene'], required=True)
    p.add_argument('--stage', choices=['hard', 'confirmation', 'merge', 'both'], default='hard')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    runner = load('original_recovery', root/'membership_stage_panel.py')
    hook = load('known_gene_noise', root/'gene_thinning_hook.py')
    rec = start_run(a.out, dict(script=Path(__file__), runner=root/'membership_stage_panel.py',
        hook=root/'gene_thinning_hook.py', design=a.task/'DESIGN.md', gate=a.task/'law_v1/report.json',
        frozen_simulator=a.task.parent/'recovery_20261107/recovery_panel_v1/manifest.json'),
        dict(mode=a.mode, stage=a.stage, oracle_gene_size=10. if a.mode == 'gene' else None),
        dict(scope='Original frozen development cases, known gene noise only; no biological calibration'), mode='gene_thinning_recovery')
    trace, previous = [], sys.argv
    try:
        assert json.loads((a.task/'law_v1/report.json').read_text())['all_checks_passed']
        defaults = json.loads((a.task.parent/'recovery_20261107/recovery_panel_v1/manifest.json').read_text())['simulator_defaults']
        assert 'phi=0.1,' in defaults
        (a.out/'folds').mkdir()
        sys.argv = [str(root/'membership_stage_panel.py'), '--out', str(a.out/'panel'), '--stage', a.stage,
            '--panel-root', str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with hook.gene_native(trace, a.out/'folds') if a.mode == 'gene' else nullcontext():
            runner.main()
        (a.out/'trace.json').write_text(json.dumps(trace, indent=2)+'\n')
    except BaseException as error:
        (a.out/'trace_partial.json').write_text(json.dumps(trace, indent=2)+'\n')
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = previous
    finish_run(a.out, rec, backend='CPU fixed full workflow; original control or known-size gene-native folds')


if __name__ == '__main__':
    main()
