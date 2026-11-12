"""Independent pooled-rate reconstruction and precision-aware HMM ablations."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.linalg import lstsq
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def lse(x,axis):
    top=np.max(x,axis=axis,keepdims=True)
    return np.squeeze(top+np.log(np.sum(np.exp(x-top),axis=axis,keepdims=True)),axis=axis)


def hmm(emissions,transition,initial,reset,dtype=np.longdouble):
    E,T,I=[np.asarray(v,dtype=dtype) for v in (emissions,transition,initial)]
    post=np.empty_like(E);path=np.empty(len(E),int)
    starts=np.flatnonzero(reset);assert starts[0]==0
    for lo,hi in zip(starts,np.r_[starts[1:],len(E)]):
        e=E[lo:hi];n=len(e);a=np.empty_like(e);b=np.zeros_like(e)
        delta=np.empty_like(e);back=np.zeros_like(e,dtype=int)
        a[0]=I+e[0];a[0]-=lse(a[0],0);delta[0]=I+e[0]
        for j in range(1,n):
            a[j]=lse(a[j-1,:,None]+T,0)+e[j];a[j]-=lse(a[j],0)
            scores=delta[j-1,:,None]+T
            back[j]=np.argmax(scores,axis=0);delta[j]=np.max(scores,axis=0)+e[j]
        for j in range(n-2,-1,-1):
            b[j]=lse(T+e[j+1,None,:]+b[j+1,None,:],1);b[j]-=max(b[j])
        joint=a+b;post[lo:hi]=np.exp(joint-lse(joint,1)[:,None])
        state=int(np.argmax(delta[-1]));path[hi-1]=state
        for j in range(n-1,0,-1):state=int(back[j,state]);path[lo+j-1]=state
    return np.asarray(post,float),path


def fixture():
    E=np.array([[.2,-.4,.7],[-.9,.3,.1],[.6,-.2,.8],[.5,.2,-.7],[.4,-.5,.6]])
    T=np.log(np.array([[.8,.1,.1],[.1,.8,.1],[.1,.1,.8]]));I=np.log([.2,.5,.3])
    reset=np.array([1,0,1,0,0],bool);scores=[];paths=list(itertools.product(range(3),repeat=5))
    for p in paths:scores.append(sum(E[j,p[j]]+(I[p[j]] if reset[j] else T[p[j-1],p[j]]) for j in range(5)))
    weights=np.exp(np.array(scores)-lse(np.array(scores),0));expected=np.zeros_like(E)
    for p,w in zip(paths,weights):
        for j,state in enumerate(p):expected[j,state]+=w
    post,path=hmm(E,T,I,reset);budget=128*np.finfo(float).eps*len(paths)*(1+np.max(np.abs(scores)))
    assert np.max(abs(post-expected))<=budget
    assert np.array_equal(path,paths[int(np.argmax(scores))])
    return dict(paths_enumerated=len(paths),max_difference=float(np.max(abs(post-expected))),budget=float(budget))


def centre(track,variance):
    weights=1/(.04+variance);order=np.argsort(track)
    cumulative=np.cumsum(weights[order],dtype=np.longdouble)
    value=track[order[np.flatnonzero(cumulative>=cumulative[-1]/2)[0]]]
    return track-value


def depth(base,gamma,L,XL):
    mean=float(XL[0]/L);variance=max(float(XL[1]/L)-mean**2,0)
    curvature=sum(float(b)*float(g)**2 for b,g in zip(base,gamma))
    return np.array([(float(g)*mean+.5*float(g)**2*variance)/(1+curvature*variance) for g in gamma])


def cores(post,path,reset,idx,bins,params):
    labels=np.array([int(s) if post[j,s]>=params['posterior'] else 1 for j,s in enumerate(path)])
    boundaries=[0]+[j for j in range(1,len(labels)) if reset[j] or labels[j]!=labels[j-1]]+[len(labels)]
    result=[]
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        if labels[lo]==1 or hi-lo<params['min_bins']:continue
        selected=bins.iloc[idx[lo:hi]]
        result.append(dict(first_bin=int(idx[lo]),last_bin=int(idx[hi-1]),chromosomes=sorted(set(selected.chromosome)),
                           minimum_support=float(np.min(post[np.arange(lo,hi),labels[lo]]))))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(exist_ok=False);trace=a.task/'trace_v1'
    assert json.loads((trace/'report.json').read_text())['all_public_outputs_identical']
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',records=trace/'comparisons.json'),{},
                  dict(scope='Conditional fixed-pool ablations,not a calibrated inference intervention'),mode='confirmation_observation_audit')
    rows=[];fixture_report=fixture();occurrences={};pool_checks=0
    gate=a.task.parent/'split_draw_stability_20261109/gate_v2'
    gate_runs=json.loads((gate/'runs.json').read_text())
    try:
        subprocess.run([sys.executable,'-I',str(Path(__file__).with_name('offset_recovery_audit.py')),
            '--out',str(a.out/'original_criteria'),'--task',str(a.task),'--arms',
            'trace_v1/replay/bin/panel','trace_v1/replay/gene/panel'],check=True)
        for item in json.loads((trace/'comparisons.json').read_text()):
            z=np.load(item['arrays']);bins=pd.read_csv(item['bins'],sep='\t')
            arm='gene' if '/replay/gene/' in item['source'] else 'bin'
            canonical=next(r for r in gate_runs if r['arm']==arm and r['override_enabled'])
            folds=[]
            for record in canonical['folds']:
                if record['case']!=item['seed']:continue
                fold=np.load(record['path'])
                if np.array_equal(fold['rows'],z['rows']):folds.append(fold)
            assert len(folds)==1;fold=folds[0]
            recorded=json.loads((trace/f'replay/{arm}/trace.json').read_text())
            model_ids=[]
            for model in recorded:
                if model['kind']!='model' or model['seed']!=item['seed']:continue
                saved=np.load(trace/f'replay/{arm}/arrays/model_{model["model"]:02d}.npz')
                if np.array_equal(saved['rows'],z['rows']):model_ids.append(model['model'])
            assert len(model_ids)==1
            groups=[r for r in recorded if r['kind']=='agglomerate' and r['model']==model_ids[0] and r['context']==item['context']]
            key=(arm,model_ids[0],item['context']);occurrence=occurrences.get(key,0);occurrences[key]=occurrence+1
            partition=np.load(groups[occurrence]['arrays']);cells,g=partition['cells'],partition['g']
            assert set(g)=={0,1}
            counts=fold['test'][cells];full=fold['full']
            expected=np.array([[sum(int(counts[i,j]) for i in range(len(g)) if int(g[i])==label)
                                for j in range(counts.shape[1])] for label in (0,1)])
            assert np.array_equal(expected,np.vstack([z['Ya'],z['Yb']]))
            exposure=full.sum(axis=1,dtype=np.float64);x=np.log(np.maximum(exposure,1));x-=x.mean()
            share=float(1-fold['share'])
            expected_L=np.array([share*sum(float(exposure[i]) for i in cells[g==label]) for label in (0,1)])
            assert np.array_equal(expected_L,np.r_[z['La'],z['Lb']])
            expected_XL=np.array([[share*sum(float(exposure[i]*x[i]**power) for i in cells[g==label])
                                   for power in (1,2)] for label in (0,1)])
            moment_budget=64*np.finfo(float).eps*len(cells)*(1+float(np.max(abs(expected_XL))))
            assert np.max(abs(expected_XL-np.vstack([z['XLa'],z['XLb']])))<=moment_budget
            pool_checks+=1
            raw=np.log((z['Ya'][0]+.5)/float(z['La'][0]))-np.log((z['Yb'][0]+.5)/float(z['Lb'][0]))
            variance=1/(z['Ya'][0]+.5)+1/(z['Yb'][0]+.5)
            shift=depth(z['base'],z['gamma'],float(z['La'][0]),z['XLa'][0])-depth(z['base'],z['gamma'],float(z['Lb'][0]),z['XLb'][0])
            programme=np.zeros_like(raw);condition=1.
            if item['programmes'] and len(z['W']):
                W=z['W'];w=1/(.04+variance);target=raw-shift;target-=np.median(target)
                penalty=1e-6*sum(float(np.dot(row*w,row)) for row in W)/len(W)
                augmented=np.vstack([np.sqrt(w)[:,None]*W.T,np.sqrt(penalty)*np.eye(len(W))])
                rhs=np.r_[np.sqrt(w)*target,np.zeros(len(W))]
                coefficient,_,_,_=lstsq(augmented,rhs,lapack_driver='gelsy')
                programme=W.T@coefficient;condition=float(np.linalg.cond(augmented))**2
            numerical=256*np.finfo(float).eps*(1+condition)*(1+float(np.abs(raw).sum()+np.abs(shift).sum()+np.abs(programme).sum()))
            assert np.max(abs(z['shift'][0]-shift-programme))<=numerical
            idx=z['idx'];reset=z['sub_reset'];d=centre((raw-shift-programme)[idx],variance[idx])
            assert np.max(abs(z['d'][0]-d))<=numerical
            assert np.max(abs(z['v'][0]-variance[idx]))<=numerical
            inside=~reset[1:];first=np.diff(d)[inside];sampling=(variance[idx][1:]+variance[idx][:-1])[inside]
            tau=item['params']['tau_floor']
            if len(first)>=10:
                mad=1.4826*np.median(abs(first-np.median(first)))
                tau=float(np.sqrt(max(tau**2,(mad**2-np.median(sampling))/2)))
            assert abs(tau-float(z['tau'][0]))<=numerical
            assert np.max(abs(z['var'][0]-tau**2-variance[idx]))<=numerical
            E=z['emit'][0];T=z['trans'];I=z['init']
            post,path=hmm(E,T,I,reset);double,_=hmm(E,T,I,reset,dtype=np.float64)
            starts=np.flatnonzero(reset);length=int(max(np.diff(np.r_[starts,len(E)])))
            magnitude=1+float(np.max(np.ptp(E,axis=1))+np.max(abs(T))+np.max(abs(I)))
            precision_budget=128*np.finfo(float).eps*length*magnitude
            assert np.max(abs(post-double))<=precision_budget
            # Normalized forward/backward three-state recursions: allowance for
            # both passes, logsumexp, normalization and exp at each chain step.
            # This is an operation/magnitude guard, not a formal libm enclosure.
            log_guard=64*np.finfo(np.float32).eps*length*magnitude
            probability_guard=float(np.expm1(2*log_guard))
            actual_error=float(np.max(abs(post-z['post'][0])))
            assert actual_error<=probability_guard
            native_cores=cores(z['post'][0],z['path'][0],reset,idx,bins,item['params'])
            assert [(r['first_bin'],r['last_bin']) for r in native_cores]==[(r['first_bin'],r['last_bin']) for r in item['results'][0]['segments']]
            row=dict(seed=item['seed'],arm=arm,context=item['context'],
                arrays=item['arrays'],tau=tau,reference_error=actual_error,probability_guard=probability_guard,
                double_extended_discrepancy=float(np.max(abs(post-double))),programme_condition=condition,
                reconstruction_guard=numerical,production_segments=item['results'][0]['segments'],ablations={})
            metadata=pd.read_csv(a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{item["seed"]}/model_cells.tsv',sep='\t')
            specimens=metadata.iloc[z['rows']].specimen.unique();assert len(specimens)==1;row['specimen']=str(specimens[0])
            target=(bins.iloc[idx].chromosome.to_numpy()=='chr12')&(bins.iloc[idx].start.to_numpy()<10000000)
            if not target.any():raise RuntimeError('Missing locked chr12 first-half territory')
            means=np.asarray(np.array([-1,0,1],np.float32)*(item['beta']*np.log(1.5)),float)
            tracks={'executed':(d,z['var'][0]),'without_programme':(centre((raw-shift)[idx],variance[idx]),z['var'][0]),
                    'without_nuisance':(centre(raw[idx],variance[idx]),z['var'][0]),'without_extra_tau':(d,variance[idx])}
            for name,(track,var) in tracks.items():
                emit=-.5*(np.log(2*np.pi*var[:,None])+(track[:,None]-means)**2/var[:,None])
                p,best=hmm(emit,T,I,reset)
                row['ablations'][name]=dict(first_half_bin_contrasts=track[target].tolist(),
                    first_half_sampling_variance=variance[idx][target].tolist(),
                    first_half_state_support=p[target].tolist(),
                    first_half_min_nonzero_support=float(np.min(p[target,0]+p[target,2])),
                    first_half_max_nonzero_support=float(np.max(p[target,0]+p[target,2])),
                    cores=cores(p,best,reset,idx,bins,item['params']),
                    numerical_boundary_overlap=bool(np.any(abs(p[target]-item['params']['posterior'])<=probability_guard)))
            rows.append(row)
        write_json(rows,a.out/'comparisons.json')
        write_json(dict(all_checks_passed=True,comparisons=len(rows),independent_pool_checks=pool_checks,original_criteria_checked=60,fixture=fixture_report,
                        production_changed=False,biological_specimens=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU independent augmented least squares and long-double HMM')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
