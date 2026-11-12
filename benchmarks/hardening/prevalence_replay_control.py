"""Fixed old predictive calculations, independently of new model fitting."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec=importlib.util.spec_from_file_location('prevalence_replay_'+name,
        Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    prior=args.task.parent/'adjacent_depth_state_20261110'
    cases=json.loads((args.task/'CASES.json').read_text())
    inputs=dict(script=Path(__file__),cases=args.task/'CASES.json',
        design=args.task/'REPLAY_DESIGN.md',predictor=Path(__file__).with_name('depth_state_predictive.py'),
        integrator=Path(__file__).with_name('state_predictive_confirmation.py'))
    for i in range(8):
        inputs[f'models{i}']=prior/f'case_{i:02d}_v1/DISCOVERY_MODELS_LOCKED.json'
        for order in (24,48,96):inputs[f'old_{i}_{order}']=prior/f'case_{i:02d}_v1/per_cell_{order}.npy'
    rec=start_run(args.out,inputs,dict(new_fit_run=False),
        dict(scope='Numerical replay diagnosis; no scientific threshold change'),mode='prevalence_replay_control')
    try:
        records=[]
        for i,case in enumerate(cases):
            D,ED,rD,g,metadata,provenance,EC,rC=load('state_predictive_confirmation').read_proposal(case)
            lock=json.loads(Path(provenance['proposal_lock']).read_text())
            strata=np.asarray(lock['parameters']['strata'],int)[metadata['selected']]
            with np.load(provenance['confirmation'])as z:C=z['counts'][metadata['selected']]
            models=json.loads(inputs[f'models{i}'].read_text())['models']
            for model in models:
                actual=load('depth_state_predictive').predict(D,C,ED,EC,rD,rC,strata,model)['per_cell']
                expected=np.load(inputs[f"old_{i}_{model['order']}"])
                records.append(dict(case=case['label'],order=model['order'],exact=bool(np.array_equal(actual,expected)),
                    changed_cells=int(np.count_nonzero(actual!=expected)),maximum_error=float(np.max(abs(actual-expected)))))
        write_json(dict(node=os.environ.get('SLURM_JOB_NODELIST'),new_fit_run=False,
            cpuinfo=Path('/proc/cpuinfo').read_text().split('\n\n')[0],
            libc_maps=[line for line in Path('/proc/self/maps').read_text().splitlines()if'libm.'in line or'libc.'in line],
            numpy_cpu_features=np._core._multiarray_umath.__cpu_features__,
            records=records,all24arrays_exact=all(r['exact']for r in records)),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU old fixed predictive calculations, no new model fit')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
