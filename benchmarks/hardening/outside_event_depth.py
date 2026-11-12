"""Research matching on counts outside a proposed genomic RNA event."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from haplodose.clones import CloneParams
from haplodose.io import write_json


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('outside_event_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def matched_weights(depth,labels,margin,n_strata=8):
    depth,labels,margin=np.asarray(depth,float),np.asarray(labels,int),np.asarray(margin,float)
    if (depth.ndim!=1 or labels.shape!=depth.shape or margin.shape!=depth.shape
            or set(labels)!={0,1} or not np.isfinite(depth).all() or not np.isfinite(margin).all()
            or np.any(margin<0) or np.any(margin>1) or n_strata<1):
        raise ValueError('Aligned binary populations, finite depths and margins in[0,1] required')
    ia,ib=np.flatnonzero(labels==0),np.flatnonzero(labels==1)
    edges=np.quantile(np.r_[depth[ia],depth[ib]],np.linspace(0,1,n_strata+1)[1:-1])
    sa,sb=np.searchsorted(edges,depth[ia]),np.searchsorted(edges,depth[ib])
    na=np.bincount(sa,weights=margin[ia],minlength=n_strata)
    nb=np.bincount(sb,weights=margin[ib],minlength=n_strata);mass=np.minimum(na,nb)
    result=np.zeros(len(depth))
    result[ia]=margin[ia]*np.divide(mass,na,out=np.zeros_like(mass),where=na>0)[sa]
    result[ib]=margin[ib]*np.divide(mass,nb,out=np.zeros_like(mass),where=nb>0)[sb]
    return result,edges


def depth_covariate(source,proposal,gene_bins,columns,mode):
    if mode not in ('full_outside','proposal_outside'):raise ValueError(mode)
    outside=~np.isin(gene_bins,columns)
    if not outside.any():raise ValueError('No genes remain outside candidate event')
    matrix=source if mode=='full_outside' else proposal
    counts=np.asarray(matrix[:,outside].sum(axis=1)).ravel()
    depth=np.log(np.maximum(counts,1.))
    return depth-depth.mean(),counts


def provider(mode):
    def weights(a,b,source,selected):
        full=np.asarray(source['source'].sum(axis=1)).ravel()
        log_full=np.log(np.maximum(full,1.));original_depth=log_full-log_full.mean()
        q=a['q'];labels=(q>.5).astype(int);margin=np.abs(2*q-1)
        old,old_edges=matched_weights(original_depth[selected],labels,margin)
        np.testing.assert_array_equal(old,b['matching_weights'])
        if mode=='uniform_full':
            depth,counts=original_depth,full
            new,edges=matched_weights(depth[selected],labels,np.ones(len(margin)))
        else:
            depth,counts=depth_covariate(source['source'],source['proposal'],source['gene_bins'],a['cols'],mode)
            new,edges=matched_weights(depth[selected],labels,margin)
        effective=[];positive=[]
        for k in (0,1):
            w=new[(labels==k)&(new>0)].astype(np.longdouble);positive.append(len(w))
            effective.append(0. if not len(w) else float(len(w)) if np.all(w==w[0]) else float(w.sum()**2/(w*w).sum()))
        minimum=CloneParams().min_cells
        return new,dict(mode=mode,original_weights_exact=True,original_edges=old_edges.tolist(),new_edges=edges.tolist(),
            original_mass=float(old.sum()),new_mass=float(new.sum()),effective_cells=effective,positive_cells=positive,
            minimum_cells=minimum,eligible=bool(min(positive)>=minimum and min(effective)>=minimum),
            zero_outside_counts=int(np.sum(counts[selected]==0)),changed_cells=int(np.sum(new!=old)),
            maximum_weight_change=float(np.max(abs(new-old))))
    return weights


def checks():
    depth=np.repeat(np.arange(4.),2);labels=np.array([0,0,0,1,1,1,1,1])
    margin=np.array([1.,.5,.5,.25,.5,.75,.25,1.])
    weights,edges=matched_weights(depth,labels,margin,2)
    exact=[Fraction(1,8),Fraction(1,16),Fraction(1,16),Fraction(1,4),Fraction(0),Fraction(0),Fraction(0),Fraction(0)]
    np.testing.assert_array_equal(weights,list(map(float,exact)));np.testing.assert_array_equal(edges,[1.5])
    uniform,_=matched_weights(depth,labels,np.ones(len(labels)),2)
    np.testing.assert_array_equal(uniform,[1/3,1/3,1/3,1.,0.,0.,0.,0.])
    assert math.fsum(uniform[labels==0])==math.fsum(uniform[labels==1])==1.
    for strata in (2,3,8):
        _,edges=matched_weights(depth,labels,margin,strata)
        ordered=sorted(map(float,depth));scalar=[]
        for k in range(1,strata):
            h=float(np.linspace(0,1,strata+1)[k])*(len(depth)-1);lo=math.floor(h);hi=math.ceil(h)
            scalar.append(ordered[lo]+(ordered[hi]-ordered[lo])*(h-lo))
        np.testing.assert_array_equal(edges,scalar)
    permutation=np.array([7,0,4,2,6,3,1,5]);other,_=matched_weights(depth[permutation],labels[permutation],margin[permutation],2)
    np.testing.assert_array_equal(other,weights[permutation])
    source=sparse.csr_matrix([[4.,5.,6.],[8.,10.,12.]])
    proposal=sparse.csr_matrix([[2.,2.,3.],[4.,5.,6.]])
    bins=np.array([0,1,2]);event=np.array([1])
    a=depth_covariate(source,proposal,bins,event,'proposal_outside')
    changed=source.copy();changed.data=changed.data*7
    b=depth_covariate(changed,proposal,bins,event,'proposal_outside')
    for before,after in zip(a,b):np.testing.assert_array_equal(before,after)
    np.testing.assert_array_equal(a[1],[5.,10.])
    np.testing.assert_array_equal(depth_covariate(source,proposal,bins,event,'full_outside')[1],[10.,20.])
    return dict(independent_fraction_stratum_mass=True,scalar_quantile_identity=True,cell_permutation=True,
        proposal_covariate_confirmation_perturbation_exact=True,all_event_genes_excluded=True,
        kernel_reference=load('averaged_count_matching').checks())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('full_outside','proposal_outside','uniform_full'),default='full_outside')
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:load('averaged_count_matching').study(args.out,args.task,numerical_checks=checks,
        weight_provider=provider(args.mode),extra_inputs=dict(outside_depth=Path(__file__),
            **(dict(margin_design=args.task/'MARGIN_DESIGN.md') if args.mode=='uniform_full' else {})))


if __name__=='__main__':main()
