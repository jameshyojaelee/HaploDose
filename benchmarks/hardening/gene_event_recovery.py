"""Whole original recovery panel with matched gene/bin proposal likelihoods."""
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
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('gene_event_'+name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('control','moment_bin','gene'),required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    (args.out/'events').mkdir()
    method, noise = load('gene_event_proposal'), load('reference_gene_noise')
    run = start_run(args.out,dict(script=Path(__file__),method=Path(method.__file__),
        noise=Path(noise.__file__),design=args.task/'DESIGN.md',
        runner=Path(__file__).with_name('membership_stage_panel.py'),
        stage=Path(__file__).with_name('membership_weighted_stages.py'),
        thinning=Path(__file__).with_name('gene_thinning_hook.py'),
        veto=Path(__file__).with_name('gene_veto_weights.py')),
        dict(mode=args.mode,depth_strata=5),dict(scope='Conditional development recovery; no error-control claim'),
        mode='gene_event_recovery')
    original_load, argv = noise.load, sys.argv
    trace = []
    def with_events(name):
        module = original_load(name)
        if name == 'gene_thinning_hook':
            original = module.gene_native
            @contextmanager
            def wrapped(*a, **k):
                with original(*a,**k) as registry:
                    with method.event_likelihood(registry,args.mode,trace,args.out/'events'):
                        yield
            module.gene_native = wrapped
        return module
    noise.load = with_events
    try:
        sys.argv = [str(noise.__file__),'--out',str(args.out/'candidate'),'--task',
            str(args.task.parent/'depth_stratified_gene_noise_20261109'),'--stage','both','--depth-strata','5']
        noise.main()
        report = json.loads((args.out/'candidate/report.json').read_text())
        if args.mode == 'control':
            prior = args.task.parent/'depth_stratified_gene_noise_20261109/both_v1'
            assert report['recovered'] == 4
            for row in report['cases']:
                before = prior/'panel'/f"seed_{row['seed']}"/'fit'
                after = args.out/'candidate/panel'/f"seed_{row['seed']}"/'fit'
                for path in before.iterdir():
                    if path.suffix == '.npy':
                        np.testing.assert_array_equal(np.load(path),np.load(after/path.name))
                    elif path.name.endswith(('.tsv','.tsv.gz')):
                        pd.testing.assert_frame_equal(pd.read_csv(path,sep='\t'),
                            pd.read_csv(after/path.name,sep='\t'),check_exact=True)
            for path in (prior/'folds').glob('*.npz'):
                with np.load(path) as old, np.load(args.out/'candidate/folds'/path.name) as new:
                    assert old.files == new.files
                    for name in old.files:
                        np.testing.assert_array_equal(old[name],new[name])
        write_json(dict(mode=args.mode,recovered=report['recovered'],cases=report['cases'],
            event_fits=len(trace),production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        write_json(trace,args.out/'event_trace.json')
        finish_run(args.out,run,backend='CPU same gene folds, alternative conditional NB event likelihoods')
    except BaseException as error:
        write_json(trace,args.out/'event_trace.json')
        finish_run(args.out,run,error=error)
        raise
    finally:
        noise.load, sys.argv = original_load, argv


if __name__ == '__main__':
    main()
