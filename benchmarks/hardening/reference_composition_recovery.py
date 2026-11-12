"""Original recovery with reference-derived gene proposal expectations."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('reference_composition_'+name,Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('control','pooled','depth'),required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    (args.out/'events').mkdir()
    (args.out/'reference').mkdir()
    method,noise,gene = load('reference_gene_composition'),load('reference_gene_noise'),load('gene_event_proposal')
    run = start_run(args.out,dict(script=Path(__file__),method=Path(method.__file__),noise=Path(noise.__file__),
        gene=Path(gene.__file__),design=args.task/'DESIGN.md',
        thinning=Path(__file__).with_name('gene_thinning_hook.py'),
        stage=Path(__file__).with_name('membership_weighted_stages.py'),
        runner=Path(__file__).with_name('membership_stage_panel.py'),
        veto=Path(__file__).with_name('gene_veto_weights.py')),
        dict(mode=args.mode,reference_regularizer_per_gene=.5),
        dict(scope='Conditional development recovery, no independent confirmation or biological validity'),
        mode='reference_composition_recovery')
    trace,reference_records,cell_records = [],[],[]
    method.estimator(noise,args.mode,reference_records,args.out/'reference')
    provider = None if args.mode == 'control' else method.region_provider(cell_records,args.out/'reference')
    original_load,argv = noise.load,sys.argv
    def with_gene(name):
        module = original_load(name)
        if name == 'gene_thinning_hook':
            original = module.gene_native
            @contextmanager
            def wrapped(*a,**k):
                with original(*a,**k) as registry:
                    with gene.event_likelihood(registry,'gene',trace,args.out/'events',region_provider=provider):
                        yield registry
            module.gene_native = wrapped
        return module
    noise.load = with_gene
    try:
        sys.argv = [str(noise.__file__),'--out',str(args.out/'candidate'),'--task',
            str(args.task.parent/'depth_stratified_gene_noise_20261109'),'--stage','both','--depth-strata','5']
        noise.main()
        report = json.loads((args.out/'candidate/report.json').read_text())
        if args.mode == 'control':
            previous = args.task.parent/'gene_event_proposals_20261109/gene_v1/candidate'
            prior_report = json.loads((previous/'report.json').read_text())
            assert report['recovered'] == 4
            assert [r['criteria'] for r in report['cases']] == [r['criteria'] for r in prior_report['cases']]
            for row in report['cases']:
                before = previous/'panel'/f"seed_{row['seed']}"/'fit'
                after = args.out/'candidate/panel'/f"seed_{row['seed']}"/'fit'
                for path in before.iterdir():
                    if path.suffix == '.npy':
                        np.testing.assert_array_equal(np.load(path),np.load(after/path.name))
                    elif path.name.endswith(('.tsv','.tsv.gz')):
                        pd.testing.assert_frame_equal(pd.read_csv(path,sep='\t'),pd.read_csv(after/path.name,sep='\t'),check_exact=True)
            for path in (previous/'folds').glob('*.npz'):
                with np.load(path) as old,np.load(args.out/'candidate/folds'/path.name) as new:
                    assert old.files == new.files
                    for name in old.files:
                        np.testing.assert_array_equal(old[name],new[name])
        write_json(dict(mode=args.mode,recovered=report['recovered'],cases=report['cases'],
            event_fits=len(trace),reference_populations=len(reference_records),
            production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        write_json(trace,args.out/'event_trace.json')
        write_json(dict(references=reference_records,cells=cell_records),args.out/'reference_records.json')
        finish_run(args.out,run,backend='CPU reference-only within-bin composition and native gene NB proposals')
    except BaseException as error:
        write_json(trace,args.out/'event_trace.json')
        write_json(dict(references=reference_records,cells=cell_records),args.out/'reference_records.json')
        finish_run(args.out,run,error=error)
        raise
    finally:
        noise.load,sys.argv = original_load,argv


if __name__ == '__main__':
    main()
