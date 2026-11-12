"""Independent quadrature/path checks and frozen-pool regional-effect experiment."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import roots_legendre
from scipy.stats import norm
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    import importlib.util
    spec=importlib.util.spec_from_file_location('regional_probe_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def integral(d,v,step,scale,sign,tolerance):
    """Numerically integrate the scalar likelihood ratio, without its closed form."""
    A=math.fsum(1/float(x) for x in v);B=math.fsum(float(x/y) for x,y in zip(d,v))
    mu=sign*step;mode=max(0.,sign*(B+mu/scale**2)/(A+1/scale**2))
    def logvalue(t):
        theta=sign*t
        likelihood=math.fsum(float((2*x*theta-theta**2)/(2*y)) for x,y in zip(d,v))
        return likelihood+float(norm.logpdf(theta,loc=mu,scale=scale))-float(norm.logcdf(sign*mu/scale))
    peak=logvalue(mode)
    def value(t):return math.exp(logvalue(t)-peak)
    one,e1=quad(value,0.,mode,epsabs=tolerance,epsrel=tolerance) if mode else (0.,0.)
    two,e2=quad(value,mode,np.inf,epsabs=tolerance,epsrel=tolerance)
    total=one+two
    return math.log(total)+peak,(e1+e2)/total


def enumerate_paths(d,v,reset,T,I,step,scale):
    """Exhaustive path sum using independently integrated run likelihoods."""
    n=len(d);evidence={};error=0.
    for lo in range(n):
        for hi in range(lo,n):
            if np.any(reset[lo+1:hi+1]):break
            for state,sign in ((0,-1),(2,1)):
                if scale==0:
                    mu=sign*step;value=math.fsum(float((2*x*mu-mu**2)/(2*y)) for x,y in zip(d[lo:hi+1],v[lo:hi+1]))
                else:
                    coarse,first=integral(d[lo:hi+1],v[lo:hi+1],step,scale,sign,1e-9)
                    value,second=integral(d[lo:hi+1],v[lo:hi+1],step,scale,sign,1e-12)
                    rounding=64*np.finfo(float).eps*(hi-lo+1)*(1+abs(coarse)+abs(value))
                    assert abs(coarse-value)<=first+second+rounding
                    error=max(error,second+rounding)
                evidence[lo,hi,state]=value
    paths=list(itertools.product(range(3),repeat=n));logs=[]
    for path in paths:
        value=math.fsum(float(I[path[j]] if reset[j] else T[path[j-1],path[j]]) for j in range(n))
        starts=[0]+[j for j in range(1,n) if reset[j] or path[j]!=path[j-1]]
        for lo,hi in zip(starts,[j-1 for j in starts[1:]]+[n-1]):
            if path[lo]!=1:value+=evidence[lo,hi,path[lo]]
        logs.append(value)
    maximum=max(logs);weights=[math.exp(x-maximum) for x in logs];total=math.fsum(weights)
    post=np.array([[math.fsum(w for p,w in zip(paths,weights) if p[j]==s)/total for s in range(3)] for j in range(n)])
    return dict(post=post,path=np.array(paths[int(np.argmax(logs))]),paths=len(paths),
                integration_log_error=error,log_evidence_ratio=math.log(total)+maximum)


def fixtures(model,reference):
    d=np.array([-.25,-.15,.3,.45,.2]);v=np.array([.01,.02,.025,.03,.01]);reset=np.array([1,0,1,0,0],bool)
    T=np.log([[.9,.07,.03],[.04,.92,.04],[.02,.06,.92]]);I=np.log([.2,.6,.2]);step=.24
    records=[]
    for scale in (0.,step):
        fitted=model.infer(d,v,reset,T,I,step,scale);expected=enumerate_paths(d,v,reset,T,I,step,scale)
        u=np.finfo(float).eps/2;gamma=(len(d)+1)*expected['paths']*u/(1-(len(d)+1)*expected['paths']*u)
        guard=2*len(d)*expected['integration_log_error']+gamma*(1+abs(expected['log_evidence_ratio']))
        assert np.max(abs(fitted['post']-expected['post']))<=guard
        assert np.array_equal(fitted['path'],expected['path'])
        if scale==0:
            emission=-.5*(d[:,None]-np.array([-step,0.,step]))**2/v[:,None]
            post,path=reference.hmm(emission,T,I,reset)
            assert np.max(abs(fitted['post']-post))<=guard and np.array_equal(fitted['path'],path)
        scaled=model.infer(3*d,9*v,reset,T,I,3*step,3*scale)
        assert np.max(abs(fitted['post']-scaled['post']))<=guard
        reversed=model.infer(-d,v,reset,T[::-1,::-1],I[::-1],step,scale)
        assert np.max(abs(fitted['post']-reversed['post'][:,::-1]))<=guard
        records.append(dict(scale=scale,paths=expected['paths'],guard=guard,
            maximum_posterior_discrepancy=float(np.max(abs(fitted['post']-expected['post'])))))
    return records


def grid_reference(d,v,T,I,step,reference):
    """Independent amplitude quadrature as an expanded ordinary HMM.

    A microscopic amplitude persists on a self-transition; a sign change
    independently draws the next run amplitude. Quadrature weights are not
    normalized, retaining the true truncated-normal density and tail bound.
    """
    n=len(d);prior_mass=float(norm.cdf(1.));null_log=float(I[1]+(n-1)*T[1,1])
    maximum_log_ratio=math.fsum(float(x*x/(2*y)) for x,y in zip(d,v))
    numerical_tail_target=1e-12
    z=max(12.,math.sqrt(2*max(0.,maximum_log_ratio-null_log+math.log(2*n/(prior_mass*numerical_tail_target)))))
    cutoff=step*(1+z)
    log_tail_bound=math.log(2*n/prior_mass)-z*z/2+maximum_log_ratio-null_log
    tail_bound=math.exp(log_tail_bound)
    records=[];previous=None;accepted=None
    for order in (32,64,128,256,512):
        nodes,weights=roots_legendre(order);amplitude=(nodes+1)*cutoff/2
        weights=weights*cutoff/2*norm.pdf(amplitude,loc=step,scale=step)/prior_mass
        assert np.isfinite(weights).all() and np.all(weights>0)
        values=np.r_[-amplitude,0.,amplitude];state=np.r_[np.zeros(order,int),1,np.full(order,2,int)]
        draw_weights=np.r_[weights,1.,weights];size=len(values)
        transition=np.full((size,size),-np.inf)
        for a in range(size):
            other=np.flatnonzero(state!=state[a])
            transition[a,other]=T[state[a],state[other]]+np.log(draw_weights[other])
            transition[a,a]=T[state[a],state[a]]
        initial=I[state]+np.log(draw_weights)
        emission=(d[:,None]*values-.5*values**2)/v[:,None]
        post,_=reference.hmm(emission,transition,initial,np.r_[True,np.zeros(n-1,bool)])
        aggregate=np.column_stack([post[:,state==s].sum(1) for s in range(3)])
        # Accumulation allowance counts both dense passes over this state space.
        u=np.finfo(np.longdouble).eps/2;operations=64*n*size;gamma=operations*u/(1-operations*u)
        magnitude=1+float(np.max(abs(emission)))+float(np.max(abs(initial)))
        input_u=np.finfo(float).eps/2;input_gamma=16*n*input_u/(1-16*n*input_u)
        output_gamma=size*input_u/(1-size*input_u)
        guard=(gamma+input_gamma)*magnitude+output_gamma+tail_bound
        difference=None if previous is None else float(np.max(abs(aggregate-previous)))
        records.append(dict(order=order,states=size,prior_quadrature_mass=float(weights.sum()),
            change_from_previous=difference,rounding_and_tail_guard=float(guard),
            rounding_and_tail_guard_extended=str(guard),tail_bound=tail_bound))
        if previous is not None and difference<=guard:
            accepted=(aggregate,guard,records)
        else:accepted=None
        previous=aggregate
    assert accepted is not None,'Amplitude quadrature has not converged at the final prespecified order'
    return accepted


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--mode',choices=('preflight','panel','audit'),required=True)
    p.add_argument('--audit-method',choices=('exhaustive','quadrature'),default='exhaustive')
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1]
    model=load('regional_amplitude_model');reference=load('confirmation_observation_audit')
    trace=root/'dev/confirmation_observation_20261109/trace_v1'
    inputs=dict(script=Path(__file__),model=Path(model.__file__),reference=Path(reference.__file__),design=a.task/'DESIGN.md')
    if a.mode!='preflight':
        inputs.update(comparisons=trace/'comparisons.json',preflight=a.task/'preflight_v1/report.json',
            previous_audit=root/'dev/confirmation_observation_20261109/audit_v2/report.json')
        assert json.loads(inputs['preflight'].read_text())['all_checks_passed']
        assert json.loads(inputs['previous_audit'].read_text())['all_checks_passed']
        items=json.loads(inputs['comparisons'].read_text())
        for i,item in enumerate(items):inputs[f'payload_{i}']=Path(item['arrays']);inputs[f'bins_{i}']=Path(item['bins'])
    if a.mode=='audit':inputs.update(panel=a.task/'panel_v1/run_manifest.json',results=a.task/'panel_v1/results.json')
    if a.audit_method=='quadrature':inputs['audit_amendment']=a.task/'AUDIT_AMENDMENT.md'
    rec=start_run(a.out,inputs,dict(mode=a.mode,scale='one existing difference step'),
        dict(scope='Fixed-pool conditional RNA contrast experiment,not whole-workflow recovery/error or DNA validation'),mode='regional_amplitude_probe')
    try:
        if a.mode=='preflight':
            checks=fixtures(model,reference)
            report=dict(all_checks_passed=True,fixtures=checks,production_changed=False)
        elif a.mode=='panel':
            rows=[]
            for i,item in enumerate(items):
                z=np.load(item['arrays']);d=z['d'][0];v=z['var'][0];step=float(np.float32(item['beta']*np.log(1.5)))
                point=model.infer(d,v,z['sub_reset'],z['trans'],z['init'],step,0.)
                means=np.array([-step,0.,step]);emit=-.5*(d[:,None]-means)**2/v[:,None]
                independent,path=reference.hmm(emit,z['trans'],z['init'],z['sub_reset'])
                n=len(d);magnitude=1+float(np.max(abs(emit)))+sum(abs(b['log_evidence_ratio']) for b in point['blocks'])
                u=np.finfo(float).eps/2;gamma=128*n*u/(1-128*n*u);guard=gamma*magnitude
                assert np.max(abs(point['post']-independent))<=guard and np.array_equal(point['path'],path)
                fitted=model.infer(d,v,z['sub_reset'],z['trans'],z['init'],step,step)
                bins=pd.read_csv(item['bins'],sep='\t');idx=z['idx']
                target=(bins.iloc[idx].chromosome.to_numpy()=='chr12')&(bins.iloc[idx].start.to_numpy()<10000000)
                assert target.any()
                np.savez_compressed(a.out/f'comparison_{i:04d}.npz',post=fitted['post'],path=fitted['path'],point_post=point['post'],point_path=point['path'])
                rows.append(dict(index=i,seed=item['seed'],context=item['context'],arm='gene' if '/replay/gene/' in item['source'] else 'bin',
                    saved_payload=item['arrays'],step=step,point_limit_reference_discrepancy=float(np.max(abs(point['post']-independent))),point_limit_guard=guard,
                    executed_first_half_max=float(np.max(z['post'][0,target,0]+z['post'][0,target,2])),
                    shared_first_half_max=float(np.max(fitted['post'][target,0]+fitted['post'][target,2])),
                    shared_first_half_min=float(np.min(fitted['post'][target,0]+fitted['post'][target,2])),
                    native_cores=reference.cores(z['post'][0],z['path'][0],z['sub_reset'],idx,bins,item['params']),
                    shared_cores=reference.cores(fitted['post'],fitted['path'],z['sub_reset'],idx,bins,item['params']),blocks=fitted['blocks']))
            (a.out/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
            report=dict(all_checks_passed=True,conditional_comparisons=len(rows),full_workflow_recovery='not run',production_changed=False)
        else:
            manifest=json.loads(inputs['panel'].read_text());assert manifest['status']=='complete'
            for group in ('input','output'):
                for record in manifest[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
            rows=json.loads(inputs['results'].read_text());assert len(rows)==len(items)==27
            for row,item in zip(rows,items):
                z=np.load(item['arrays']);pred=np.load(a.task/f'panel_v1/comparison_{row["index"]:04d}.npz')
                bins=pd.read_csv(item['bins'],sep='\t')
                expected=reference.cores(pred['post'],pred['path'],z['sub_reset'],z['idx'],bins,item['params'])
                assert expected==row['shared_cores']
            # Prespecified first production109 first-half record, before outcomes.
            chosen=next(r for r in rows if r['seed']==20261109 and r['arm']=='bin')
            assert chosen['index']==9 and chosen['context'].endswith('bins209-217')
            item=items[chosen['index']];z=np.load(item['arrays']);bins=pd.read_csv(item['bins'],sep='\t')
            idx=np.flatnonzero(bins.iloc[z['idx']].chromosome.to_numpy()=='chr12');assert len(idx)>0
            assert np.array_equal(idx,np.arange(idx[0],idx[-1]+1))
            d,v=z['d'][0,idx],z['var'][0,idx];reset=z['sub_reset'][idx].copy();assert reset[0]
            pred=np.load(a.task/f'panel_v1/comparison_{chosen["index"]:04d}.npz');actual=pred['post'][idx]
            if a.audit_method=='exhaustive':
                assert len(idx)<=12
                expected=enumerate_paths(d,v,reset,z['trans'],z['init'],chosen['step'],chosen['step'])
                n=len(idx);u=np.finfo(float).eps/2;operations=(n+1)*expected['paths'];gamma=operations*u/(1-operations*u)
                guard=2*n*expected['integration_log_error']+gamma*(1+abs(expected['log_evidence_ratio']))
                expected_post=expected['post'];details=dict(exhaustive_actual_chr12_paths=expected['paths'])
                assert np.array_equal(pred['path'][idx],expected['path'])
            else:
                assert not np.any(reset[1:]),'Audit must retain one whole original reset block'
                expected_post,guard,orders=grid_reference(d,v,z['trans'],z['init'],chosen['step'],reference)
                details=dict(actual_chr12_positions=len(idx),quadrature_orders=orders,
                    path_validation='MAP path checked by exhaustive preflight and all27 point-limit references; actual full-block audit targets posterior')
            assert np.max(abs(actual-expected_post))<=guard
            report=dict(all_checks_passed=True,conditional_comparisons_recounted=27,selected_comparison=chosen['index'],
                quadrature_posterior_discrepancy=float(np.max(abs(actual-expected_post))),
                quadrature_posterior_guard=float(guard),quadrature_posterior_guard_extended=str(guard),
                **details,production_changed=False,biological_specimens=0)
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        finish_run(a.out,rec,backend='CPU analytic integration and independent quadrature/exhaustive paths')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
