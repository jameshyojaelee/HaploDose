"""Original five recovery cases with conditional-objective diagnostics."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
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
    p.add_argument('--arm', choices=['legacy', 'conditional', 'NB'], required=True)
    p.add_argument('--trace-applied', action='store_true')
    p.add_argument('--offset-aware-profiles', action='store_true')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    runner = load('recovery', root/'membership_stage_panel.py')
    helper = load('nb_offset', root/'offset_nb_objective.py')
    profile_hook = load('profile_offset', root/'profile_offset.py') if a.offset_aware_profiles else None
    inputs = dict(script=Path(__file__), runner=root/'membership_stage_panel.py',
        helper=root/'offset_nb_objective.py', design=a.task/'DESIGN.md')
    params = dict(arm=a.arm)
    if profile_hook:
        inputs['profile_hook'] = root/'profile_offset.py'
        params['offset_aware_profiles'] = True
    if a.trace_applied:
        params['trace_applied'] = True
    rec = start_run(a.out, inputs, params,
        dict(scope='Conditional-objective diagnostic, original five development cases; no biological calibration'),
        mode='offset_objective_panel')
    trace, worst, previous = [], {}, sys.argv
    profile_trace, profile_contexts = [], {}
    try:
        sys.argv = [str(root/'membership_stage_panel.py'), '--out', str(a.out/'panel'), '--stage', 'hard',
            '--panel-root', str(a.task.parent/'recovery_20261107/recovery_panel_v1'),
            '--assertion-source', str(a.task.parent/'recovery_20261107/baseline/tests/test_split_mode.py')]
        with helper.objective_trace(a.arm, trace, worst, trace_applied=a.trace_applied), \
                (profile_hook.conditioned_profiles(profile_trace, profile_contexts) if profile_hook else nullcontext()):
            runner.main()
        if profile_hook:
            (a.out/'profile_trace.json').write_text(json.dumps(profile_trace, indent=2)+'\n')
            expected_contexts = []
            for (model, nonzero), values in sorted(profile_contexts.items()):
                path = a.out/f'expected_model_{model}_nonzero_{int(nonzero)}.npz'
                np.savez_compressed(path, **values)
                expected_contexts.append(dict(model=model, nonconstant_offset=nonzero, path=str(path)))
            (a.out/'expected_contexts.json').write_text(json.dumps(expected_contexts, indent=2)+'\n')
        (a.out/'trace.json').write_text(json.dumps(trace, indent=2)+'\n')
        contexts = []
        for (model, kind), item in sorted(worst.items()):
            path = a.out/f'worst_model_{model}_{kind}.npz'
            np.savez_compressed(path, **item['data'], candidate_offset=item['candidate_offset'],
                candidate_doses=item['candidate_doses'], **item.get('extra', {}))
            contexts.append(dict(path=str(path), row=item['row']))
        (a.out/'contexts.json').write_text(json.dumps(contexts, indent=2)+'\n')
        summary = []
        for kind in (('profile', 'profile_applied', 'offset') if a.trace_applied else ('profile', 'offset')):
            rows = [v for v in trace if v['kind'] == kind]
            summary.append(dict(kind=kind, updates=len(rows), decreased=sum(v['decreased'] for v in rows),
                min_gain=min(v['gain'] for v in rows),
                max_gain=max(v['gain'] for v in rows),
                optimizer_failures=sum(not v['optimizer']['success'] for v in rows if 'optimizer' in v)))
        (a.out/'report.json').write_text(json.dumps(dict(arm=a.arm, summary=summary,
            worst_contexts=len(contexts), biological_specimens=0, production_changed=False,
            scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finally:
        sys.argv = previous
    finish_run(a.out, rec, backend='CPU original recovery and float64 conditional NB diagnostics')


if __name__ == '__main__':
    main()
