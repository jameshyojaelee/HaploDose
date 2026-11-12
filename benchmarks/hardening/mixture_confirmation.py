"""Experimental observable-mixture emissions; no package default or calibration claim."""
from contextlib import contextmanager
import importlib.util
from pathlib import Path
import sys
import numpy as np
from scipy import sparse
from haplodose import clones, profile


def allocation_moments(probability, exposure, weights):
    q=np.asarray(probability,float);L=np.asarray(exposure,float)
    if q.shape!=L.shape or q.ndim!=1 or not np.isfinite(q).all() or np.any((q<0)|(q>1)):
        raise ValueError('Invalid aligned probabilities/exposures')
    if not np.isfinite(L).all() or np.any(L<0):raise ValueError('Invalid exposure')
    W=sparse.csr_matrix(weights,dtype=float)
    if W.shape[0]!=len(q) or np.any(W.data<0) or not np.isfinite(W.data).all():raise ValueError('Invalid weights')
    A=W.multiply(L[:,None]);total=np.asarray(A.sum(axis=0)).ravel()
    if np.any(total<=0):raise ValueError('Zero exposure')
    p=np.asarray(A.T@q).ravel()/total
    h=np.asarray(A.multiply(A).T@(q*(1-q))).ravel()/total**2
    return p,h


def observable_states(pa,pb,ha,hb,delta):
    pa,pb,ha,hb=[np.atleast_1d(np.asarray(v,float))[:,None] for v in (pa,pb,ha,hb)]
    d=np.asarray(delta,float)[None]
    e=np.expm1(d);a=1+pa*e;b=1+pb*e
    mean=np.log1p(pa*e)-np.log1p(pb*e)
    mean=np.where((pa==1)&(pb==0),d,mean)
    mean=np.where((pa==0)&(pb==1),-d,mean)
    variance=ha*(e/a)**2+hb*(e/b)**2
    return mean,variance


def load_control():
    path=Path(__file__).with_name('membership_weighted_stages.py')
    spec=importlib.util.spec_from_file_location('merge_weight_control',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


@contextmanager
def weighted_stages(stage,trace):
    if stage not in {'merge_control','mixture'}:raise ValueError(stage)
    with load_control().weighted_stages('merge',trace):
        if stage=='merge_control':
            yield
            return
        if sys.getprofile() is not None:raise RuntimeError('Existing profiler')
        M=clones.SpecimenModel
        event,score,differences=M.event_mixtures,M._score_candidate,M.differences
        # The control wraps event_mixtures, so observe its original function's
        # locals at the _fit_event_mixture return instead of matching wrapper code.
        cache={};active=None
        def events(self,*args,**kwargs):
            if self.alleles is not None or self.graph is not None:raise ValueError('Expression-only nonspatial prototype')
            cells=kwargs.get('cells',args[5] if len(args)>5 else None);captured={}
            def observe(frame,kind,value):
                if kind!='return' or frame.f_code is not clones._fit_event_mixture.__code__ or not isinstance(value,dict):return
                loc=frame.f_back.f_locals
                if loc.get('self') is not self or 'cols' not in loc or 'part' not in loc:return
                cols,part=loc['cols'],loc['part']
                source=f"event_chrom{loc['c']}_bins{cols[part[0]]}-{cols[part[-1]]}"
                captured[source]=value['responsibility'].copy()
            sys.setprofile(observe)
            try:found=event(self,*args,**kwargs)
            finally:sys.setprofile(None)
            if cells is None:raise ValueError('Missing cell identities')
            key=(id(self),np.asarray(cells).tobytes());cache[key]={}
            for g,source in found:
                q=captured[source]
                if not np.array_equal(q>.5,g.astype(bool)):raise ValueError('Probability identity mismatch')
                cache[key][source]=q
            return found
        def scored(self,k,cells,g,source):
            nonlocal active
            if not source.startswith('event_'):return score(self,k,cells,g,source)
            q=cache[(id(self),np.asarray(cells).tobytes())][source]
            if g.max()!=1:raise ValueError('Binary proposal required')
            fullq=np.zeros(self.rows.n);fullq[cells]=1-q
            W=sparse.csr_matrix((np.ones(len(cells)),(cells,g)),shape=(self.rows.n,2))
            before=active;active=(self,fullq,W,source)
            entry=dict(stage='mixture_confirmation',source=source,rows=self.rows.rows.tolist(),calls=[])
            active=active+(entry,)
            try:result=score(self,k,cells,g,source)
            finally:active=before
            entry['accepted']=result is not None;trace.append(entry)
            return result
        def difference(self,*args,**kwargs):
            if active is None or active[0] is not self:return differences(self,*args,**kwargs)
            _,q,W,source,entry=active
            caller=sys._getframe(1)
            if caller.f_code is M.agglomerate.__code__:
                if caller.f_locals['A']!=[0] or caller.f_locals['B']!=[1]:raise ValueError('Unexpected agglomeration pair')
                context='initial'
            elif caller.f_code.co_name=='matched_differences' and 'W' in caller.f_locals:
                W=caller.f_locals['W'];context='depth_matched'
            else:raise ValueError('Unrecognized confirmation pooling path')
            p,h=allocation_moments(q,self.L,W)
            if len(p)!=2:raise ValueError('Expected one binary pair')
            if p[0]==p[1]:
                return [dict(segments=[],distance=float('inf'),testable=False,untestable_reason='equal_mixture_purities')]
            gaussian=profile.gaussian_logpdf
            def emit(x,mean,var):
                states=np.asarray(mean)
                if states.shape!=(3,):raise ValueError('Unexpected difference state grid')
                mixed,extra=observable_states(p[:1],p[1:],h[:1],h[1:],states)
                return gaussian(x,mixed[:,None,:],var+extra[:,None,:])
            profile.gaussian_logpdf=emit
            try:result=differences(self,*args,**kwargs)
            finally:profile.gaussian_logpdf=gaussian
            entry['calls'].append(dict(context=context,purity=p.tolist(),allocation_factor=h.tolist(),
                                       segments=result[0]['segments']))
            return result
        M.event_mixtures=events;M._score_candidate=scored;M.differences=difference
        try:yield
        finally:
            M.event_mixtures=event;M._score_candidate=score;M.differences=differences
