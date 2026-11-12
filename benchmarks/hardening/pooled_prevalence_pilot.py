"""Fresh complete-workflow test of the frozen shared-prevalence safeguard."""
import argparse
from contextlib import nullcontext
from functools import lru_cache
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


PHASE=8
ARMS=('production','pooled_prevalence')


def frozen_root():
    return Path(os.environ['HAPLODOSE_DATA_ROOT'])/'dev/pooled_prevalence_workflow_20261110/source_v2'


@lru_cache(None)
def load(name):
    path=frozen_root()/'benchmarks/hardening'/(name+'.py')
    spec=importlib.util.spec_from_file_location('frozen_pooled_pilot_'+name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def verify_frozen():
    root=frozen_root()
    manifest=json.loads((root/'source_manifest.json').read_text())
    for relative,digest in manifest.items():
        if hashlib.sha256((root/relative).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Changed frozen full-workflow source: '+relative)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--condition',choices=('matched','state','two_clone'),required=True)
    parser.add_argument('--start',type=int,required=True)
    parser.add_argument('--stop',type=int,required=True)
    args=parser.parse_args();assert 0<=args.start<args.stop<=10
    args.out.mkdir(exist_ok=False)
    workflow=args.task.parent/'pooled_prevalence_workflow_20261110'
    inputs=dict(script=Path(__file__),design=args.task/'PILOT_DESIGN.md',
        full_gate=workflow/'candidate_v2/report.json',
        full_manifest=workflow/'candidate_v2/run_manifest.json',
        strict_control=workflow/'strict_control_v1/report.json',
        frozen_manifest=frozen_root()/'source_manifest.json')
    for name in('gene_veto_fresh','adjacent_reference_noise','gene_thinning_hook',
            'gene_shared_confirmation','gene_veto_weights','pooled_prevalence_workflow',
            'profiled_prevalence','depth_state_predictive','state_predictive_confirmation'):
        inputs[name]=frozen_root()/'benchmarks/hardening'/(name+'.py')
    # Gate BEFORE any generation, even when an upstream process exits0 while
    # reporting a failed original assertion.
    report=json.loads(inputs['full_gate'].read_text())
    assert json.loads(inputs['full_manifest'].read_text())['status']=='complete'
    assert report['recovered']==5 and report['full_workflow_including_merge']
    assert all(r['original_assertion_passed']and all(c['passed']for c in r['criteria'].values())
               for r in report['cases'])
    assert json.loads(inputs['strict_control'].read_text())['strict_control']['exact']
    verify_frozen()
    rec=start_run(args.out,inputs,dict(phase=PHASE,condition=args.condition,
        start=args.start,stop=args.stop,arms=ARMS),
        dict(scope='Fresh conditional complete-workflow simulation; biological n=0'),
        mode='pooled_prevalence_pilot')
    fresh,noise,hook,method,veto,guard=[load(name)for name in('gene_veto_fresh',
        'adjacent_reference_noise','gene_thinning_hook','gene_shared_confirmation',
        'gene_veto_weights','pooled_prevalence_workflow')]
    params=clones.CloneParams(split_mode='thin')
    results=[]
    try:
        for rep in range(args.start,args.stop):
            seed,X,genes,metadata,positions,arrays=fresh.generate(rep,PHASE,args.condition)
            case=args.out/f'seed_{seed}';case.mkdir()
            fresh.write_input(case/'input',X,genes,metadata,positions)
            np.savez_compressed(case/'generation.npz',**arrays)
            for arm in ARMS:
                out=case/arm;out.mkdir()
                for name in('folds','weights','arrays','guards'):(out/name).mkdir()
                estimates,split_trace,stage_trace,veto_trace,guard_trace=[],[],[],[],[]
                candidate=arm=='pooled_prevalence'
                with veto.observe_weights(veto_trace,out/'weights',apply=True,weight_units='effective')if candidate else nullcontext():
                    with hook.gene_native(split_trace,out/'folds',gene_size=noise.estimator(estimates,5),
                            case_seeds=(seed,))if candidate else nullcontext()as registry:
                        with guard.safeguard(method,registry,'depth_gene_score',stage_trace,out/'arrays',
                                guard_trace,out/'guards')if candidate else nullcontext():
                            pipeline.run(case/'input/counts.h5',case/'input/model_cells.tsv',
                                case/'input/positions.tsv',out/'fit',params=params,chunk=97,
                                counts_kind='raw',genome_build='hg38',
                                coordinate_system='one-based-closed',gene_namespace='symbol')
                write_json(dict(estimates=estimates,splits=split_trace,stages=stage_trace,
                    veto=veto_trace,guards=guard_trace),out/'trace.json')
                row=fresh.score(out/'fit',metadata,args.condition)
                row.update(replicate=rep,phase=PHASE,seed=seed,arm=arm,case=str(case),
                    fit_manifest=file_record(out/'fit/run_manifest.json'))
                results.append(row);write_json(results,args.out/'results.json')
                print(json.dumps(row),flush=True)
        verify_frozen()
        finish_run(args.out,rec,backend='CPU paired production/frozen shared-prevalence whole workflows')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
