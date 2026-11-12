"""Locate processor-dependent proposal decisions using identical saved inputs."""
import argparse
from contextlib import contextmanager,nullcontext
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numpy as np
from haplodose import clones
from haplodose.provenance import start_run,finish_run,file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


@contextmanager
def traced_search(helper,directory):
    entries=[]
    with helper.capture_search() as predictive:
        original=sys.getprofile()
        def observe(frame,event,result):
            original(frame,event,result)
            if event!='return':return
            if frame.f_code is clones._fit_event_mixture.__code__:
                outer=frame.f_back
                if outer.f_code is not clones.SpecimenModel.event_mixtures.__code__:return
                loc=outer.f_locals;m=loc['self'];cols=loc['cols'][loc['part']]
                data={key:np.asarray(frame.f_locals[key]).copy() for key in ('O','E','r')}
                data.update(cells=np.asarray(loc['cells']).copy(),cols=cols.copy())
                record=dict(stage='event_fit',index=len(entries),first_bin=int(cols[0]),last_bin=int(cols[-1]),
                            returned=result is not None)
                fit=frame.f_locals.get('fit')
                if fit is not None:
                    record.update(solver_success=bool(fit.success),solver_message=str(fit.message))
                    data['optimizer_x']=fit.x.copy()
                if result is not None:
                    data['q']=result['responsibility'].copy()
                    q=data['q'];sizes=[int(np.sum(q<=.5)),int(np.sum(q>.5))]
                    record.update(gain=float(result['gain']),budget=float(result['numerical_gain_budget']),sizes=sizes,
                        eligible=bool(result['gain']>result['numerical_gain_budget'] and min(sizes)>=m.params.min_cells))
                else:
                    for key in ('gradient','projected','score_budget','curvature','components'):
                        if key in frame.f_locals:data[key]=np.asarray(frame.f_locals[key]).copy()
                    for key in ('budget','local_gain_error','curvature_ok'):
                        if key in frame.f_locals:record[key]=float(frame.f_locals[key])
            elif frame.f_code is clones.proposal_features.__code__:
                if result is None:return
                data=dict(F=result.copy())
                record=dict(stage='proposal_features',index=len(entries))
            elif frame.f_code is clones.leiden_groups.__code__:
                if result is None:return
                loc=frame.f_locals
                data=dict(F=loc['F'].copy(),groups=np.asarray(result).copy())
                for name in ('X','idx','edges'):
                    if name in loc:data[name]=loc[name].copy()
                if 'part' in loc:data['raw_membership']=np.asarray(loc['part'].membership)
                record=dict(stage='leiden_groups',index=len(entries))
            elif frame.f_code is clones.SpecimenModel.em.__code__:
                if result is None:return
                data=dict(R=result[2].copy(),doses=np.asarray(result[0]),pi=np.asarray(result[1]))
                record=dict(stage='em',index=len(entries))
            else:return
            path=directory/f'{len(entries):04d}_{record["stage"]}.npz'
            np.savez_compressed(path,**data);record['arrays']=file_record(path);entries.append(record)
        sys.setprofile(observe)
        try:yield entries,predictive
        finally:sys.setprofile(original)


def comparison(args):
    report=[]
    for replicate in (150,167):
        paths=[args.task/f'{host}_{args.attempt}/replicate_{replicate}' for host in ('intel','amd')]
        summaries=[json.loads((p/'summary.json').read_text()) for p in paths]
        for host in ('intel','amd'):
            manifest=json.loads((args.task/f'{host}_{args.attempt}/run_manifest.json').read_text())
            if manifest['status']!='complete':raise ValueError('Incomplete processor replay')
        with np.load(paths[0]/'result.npz') as x,np.load(paths[1]/'result.npz') as y:
            identity=np.array_equal(x['R'],y['R'])
            maximum=None if x['R'].shape!=y['R'].shape else float(np.max(np.abs(x['R']-y['R'])))
        traces=[s['trace'] for s in summaries];first=None;changed=[]
        for i,(a,b) in enumerate(zip(*traces)):
            for row in (a,b):
                if file_record(row['arrays']['path'])['sha256']!=row['arrays']['sha256']:raise ValueError('Trace changed')
            difference={}
            for key in set(a)|set(b):
                if key!='arrays' and a.get(key)!=b.get(key):difference[key]=[a.get(key),b.get(key)]
            with np.load(a['arrays']['path']) as x,np.load(b['arrays']['path']) as y:
                for name in set(x.files)|set(y.files):
                    if name not in x or name not in y:difference[name]='missing on one processor'
                    elif x[name].shape!=y[name].shape:difference[name]=dict(shapes=[list(x[name].shape),list(y[name].shape)])
                    elif not np.array_equal(x[name],y[name]):
                        difference[name]=dict(maximum_absolute_difference=float(np.max(np.abs(x[name]-y[name]))))
            if difference:
                entry=dict(index=i,stage=a['stage'],differences=difference)
                if first is None:first=entry
                if any(k in difference for k in ('eligible','returned','first_bin','last_bin')):changed.append(entry)
        report.append(dict(replicate=replicate,identical_final_probabilities=identity,maximum_final_difference=maximum,
            trace_lengths=[len(t) for t in traces],first_divergence=first,changed_decisions=changed,
            predictive=[s['predictive'] for s in summaries],old_same_processor_identity=[s['historical_R_exact'] for s in summaries]))
    args.out.mkdir(exist_ok=False)
    (args.out/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=['intel','amd','compare'],required=True)
    parser.add_argument('--previous',type=Path,required=True)
    parser.add_argument('--attempt',default='v1')
    parser.add_argument('--precision-hook',action='store_true')
    args=parser.parse_args()
    if args.mode=='compare':comparison(args);return
    args.out.mkdir(exist_ok=False)
    old_source=args.previous/'predictive_rank_source_v3/benchmarks/hardening'
    rank=load('frozen_rank',old_source/'global_rank_gate.py')
    generator=rank.load_generator(old_source/'global_rank_generator.py')
    helper_path=Path(__file__).with_name('membership_predictive_search.py')
    helper=load('predictive_search',helper_path)
    historic=args.previous/('predictive_rank_v4' if args.mode=='intel' else 'predictive_rank_v3')/'evaluation.jsonl'
    previous={r['replicate']:r for r in map(json.loads,historic.read_text().splitlines())
              if r['condition']=='two_clones' and r['replicate'] in (150,167)}
    params=clones.CloneParams(n_programmes=0,spatial_neighbors=0,seed=20261107)
    inputs=dict(runner=__file__,observer=helper_path,previous=historic,
                rank=old_source/'global_rank_gate.py',generator=old_source/'global_rank_generator.py')
    correction=None
    if args.precision_hook:
        inputs['precision_hook']=Path(__file__).with_name('precise_proposal_graph.py')
        correction=load('precision_hook',inputs['precision_hook'])
    if args.mode=='amd':
        for replicate in (150,167):inputs[str(replicate)]=args.task/f'intel_{args.attempt}/replicate_{replicate}/input.npz'
    rec=start_run(args.out,inputs,asdict(params),dict(scope='Two frozen processor-sensitive development specimens',
        mode=args.mode,host=platform.node(),precision_hook=args.precision_hook,
        input_policy='Generate once on Intel; replay exact saved arrays on AMD'),mode='hardware_trace')
    try:
        reports=[]
        for replicate in (150,167):
            directory=args.out/f'replicate_{replicate}';directory.mkdir()
            if args.mode=='intel':
                rng=np.random.default_rng(np.random.SeedSequence([20261107,20,2,replicate]))
                Y,base,phi,truth,latent=generator.generate(rng,'two_clones',600,240,return_latent=True)
                reset=np.zeros(240,bool);reset[::60]=True
                data=dict(counts=Y,reference=base,phi=phi,truth=truth,reset=reset,
                          state_effect=latent['state_effect'],mean=latent['mean'])
            else:
                with np.load(inputs[str(replicate)]) as z:data={name:z[name].copy() for name in z.files}
                Y,base,phi,reset=(data[name] for name in ('counts','reference','phi','reset'))
            hashes={name:rank.array_hash(value) for name,value in data.items()}
            # Latent mean generation differed historically between CPUs. It is
            # not passed to inference; now even that array is shared exactly.
            for name,digest in hashes.items():
                if name!='mean' and digest!=previous[replicate]['input_sha256'][name]:raise ValueError('Historical inference input differs')
            np.savez_compressed(directory/'input.npz',**data)
            def fit():
                return clones.fit_specimen(clones.Rows(Y,np.arange(len(Y))),Y.sum(axis=1,dtype=np.float64),
                                           base[:,None],reset,phi,params)
            with nullcontext() if correction is None else correction.precise_graph():
                if correction is not None:correction.checks()
                unobserved=fit()
                trace_dir=directory/'trace';trace_dir.mkdir()
                with traced_search(helper,trace_dir) as (trace,predictive):result=fit()
            np.testing.assert_array_equal(unobserved['R'],result['R'])
            np.testing.assert_array_equal(unobserved['doses'],result['doses'])
            old=previous[replicate]['assignments']
            if file_record(old['path'])['sha256']!=old['sha256']:raise ValueError('Historical assignments changed')
            with np.load(old['path']) as z:historical=np.array_equal(z['raw_R'],result['R'])
            np.savez_compressed(directory/'result.npz',R=result['R'],doses=result['doses'])
            row=dict(replicate=replicate,input_sha256=hashes,observer_identity=True,historical_R_exact=historical,
                     trace=trace,predictive=predictive,model_log=result['model'].log)
            (directory/'summary.json').write_text(json.dumps(row,indent=2,default=lambda x:x.item()))
            reports.append(dict(replicate=replicate,observer_identity=True,historical_R_exact=historical,
                                events=len(trace),eligible_proposals=len(predictive)))
            print(json.dumps(reports[-1]),flush=True)
        (args.out/'summary.json').write_text(json.dumps(reports,indent=2))
        finish_run(args.out,rec,backend='cpu')
    except Exception as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
