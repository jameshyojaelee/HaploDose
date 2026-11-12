"""Conditional exact short-event confirmation ablations, never package defaults.

The NB model remains conditional on fitted discovery memberships/nuisances.
Exact integration is not independent confirmation or biological calibration.
Long events retain the existing procedure, with that fallback explicitly logged.
"""
from contextlib import contextmanager
import hashlib
import importlib.util
from pathlib import Path
import sys
import numpy as np
from haplodose import clones, profile


def load(name):
    path=Path(__file__).with_name(name+'.py')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def balanced_subset(depth, labels, identities, n_strata=8):
    """Equal integer counts in the original depth strata, selected without RNA outcomes.

    Hash order is a prespecified reproducible pseudo-random cell ordering. It
    cannot depend on confirmation counts, q magnitude, or biological truth.
    """
    depth,labels,identities=map(np.asarray,(depth,labels,identities))
    if depth.shape!=labels.shape or depth.shape!=identities.shape or len(np.unique(identities))!=len(depth):
        raise ValueError('Invalid cell alignment')
    edges=np.quantile(depth,np.linspace(0,1,n_strata+1)[1:-1])
    strata=np.searchsorted(edges,depth);selected=[];sizes=[]
    order=[hashlib.sha256(('HaploDose-exact-depth-v1:'+str(row)).encode()).digest() for row in identities]
    for s in range(n_strata):
        groups=[np.flatnonzero((strata==s)&(labels==k)) for k in (0,1)]
        n=min(map(len,groups));sizes.append(n)
        for group in groups:selected.extend(sorted(group,key=lambda i:order[i])[:n])
    return np.sort(np.asarray(selected,int)),sizes


@contextmanager
def weighted_stages(stage,trace,*,allow_unequal=False):
    if stage not in {'merge_control','exact_initial','exact_matched'}:raise ValueError(stage)
    model=load('shared_membership_hmm')
    with load('membership_weighted_stages').weighted_stages('merge',trace):
        if stage=='merge_control':
            yield
            return
        if sys.getprofile() is not None:raise RuntimeError('Existing profiler')
        M=clones.SpecimenModel;event,score,differences=M.event_mixtures,M._score_candidate,M.differences
        cache={};active=None
        def events(self,*args,**kwargs):
            if self.alleles is not None or self.graph is not None or self.params.split_mode!='thin':
                raise ValueError('Expression-only nonspatial thinning prototype')
            if self.fold_scale['prop']!=self.fold_scale['test'] and not allow_unequal:
                raise ValueError('Unequal thinning requires explicit experimental opt-in')
            cells=kwargs.get('cells',args[5] if len(args)>5 else None);captured={}
            def observe(frame,kind,value):
                if kind!='return' or frame.f_code is not clones._fit_event_mixture.__code__ or not isinstance(value,dict):return
                loc=frame.f_back.f_locals
                if loc.get('self') is not self or 'cols' not in loc or 'part' not in loc:return
                cols=loc['cols'][loc['part']]
                source=f"event_chrom{loc['c']}_bins{cols[0]}-{cols[-1]}"
                captured[source]=dict(q=value['responsibility'].copy(),E=loc['E'].copy(),r=loc['r_nb'].copy(),cols=cols.copy())
            sys.setprofile(observe)
            try:found=event(self,*args,**kwargs)
            finally:sys.setprofile(None)
            if cells is None:raise ValueError('Missing cell identities')
            key=(id(self),np.asarray(cells).tobytes());cache[key]={}
            for g,source in found:
                item=captured[source]
                if not np.array_equal(item['q']>.5,g.astype(bool)):raise ValueError('Probability identity mismatch')
                cache[key][source]=item
            return found
        def scored(self,k,cells,g,source):
            nonlocal active
            if not source.startswith('event_'):return score(self,k,cells,g,source)
            item=cache[(id(self),np.asarray(cells).tobytes())][source]
            entry=dict(stage='exact_confirmation',mode=stage,source=source,parent=int(k),
                       specimen_rows=self.rows.rows.tolist(),proposal_cells=np.asarray(cells).tolist(),calls=[])
            trace.append(entry)
            if 3**len(item['cols'])>1000000:
                entry.update(fallback='declared_exact_state_cap',accepted=None)
                result=score(self,k,cells,g,source);entry['accepted']=result is not None
                return result
            if g.max()!=1:raise ValueError('Binary proposal required')
            before=active;active=(self,item,np.asarray(cells),g,entry)
            try:result=score(self,k,cells,g,source)
            finally:active=before
            entry['accepted']=result is not None
            return result
        def difference(self,*args,**kwargs):
            if active is None or active[0] is not self:return differences(self,*args,**kwargs)
            _,item,cells,g,entry=active;caller=sys._getframe(1)
            if caller.f_code is M.agglomerate.__code__:
                if caller.f_locals['A']!=[0] or caller.f_locals['B']!=[1]:raise ValueError('Unexpected agglomeration pair')
                context='initial'
            elif caller.f_code.co_name=='matched_differences' and 'W' in caller.f_locals:
                context='depth_matched'
            else:raise ValueError('Unrecognized confirmation path')
            # Preserve the old ranking statistic/variance units; replace only
            # the posterior-supported intervals. It is not a new calibrated z.
            saved={}
            def observe(frame,kind,value):
                if kind=='return' and frame.f_code is profile.difference_segments.__code__:
                    saved.update({k:frame.f_locals[k].copy() for k in ('d','var','idx')})
            sys.setprofile(observe)
            try:old=differences(self,*args,**kwargs)
            finally:sys.setprofile(None)
            if len(old)!=1 or saved['d'].shape[0]!=1:raise ValueError('Single binary contrast required')
            cols=item['cols'];q,E,r=item['q'],item['E'],item['r']
            # Captured moments describe the proposal share. The confirmation
            # NB has its own share of both mean and convolution shape.
            scale=self.fold_scale['test']/self.fold_scale['prop']
            if scale!=1:
                E=E*scale;r=r*scale
            if not np.all(np.diff(cols)==1) or np.any(self.reset[cols[1:]]):raise ValueError('Event crosses discontinuity')
            choose=np.arange(len(cells));sizes=None
            if context=='depth_matched':
                choose,sizes=balanced_subset(self.x[cells],g,self.rows.rows[cells])
            call=dict(context=context,original_segments=old[0]['segments'],selected_cells=choose.tolist(),
                      per_stratum_per_group=sizes,confirmation_to_proposal_scale=float(scale))
            entry['calls'].append(call)
            if len(choose)<2*self.params.min_cells:
                call.update(exact_segments=[],reason='insufficient_depth_overlap')
                if context=='depth_matched' and stage=='exact_initial':return old
                return [dict(old[0],segments=[],testable=False,untestable_reason='insufficient_depth_overlap')]
            Y=self.rows.take(cells,self.test_fold)[:,cols]
            pi=np.clip((E.T@q)/E.sum(axis=0),0.,1.)
            fitted=model.build(Y[choose],E[choose],r[choose],q[choose],profile.DIFF_LEVELS*self.beta*np.log(1.5),
                profile.DIFF_PRIOR,self.diff_pp.p_change,left_neutral=bool(not self.reset[cols[0]]),
                right_neutral=bool(cols[-1]+1<self.J and not self.reset[cols[-1]+1]),mixture_proportion=pi)
            fit=model.enumerate_states(fitted)
            reset=np.r_[True,np.zeros(len(cols)-1,bool)];segments=[]
            for a,b,s in profile.credible_difference_runs(fit['path'],fit['posterior'],reset,self.diff_pp):
                first,last=int(cols[a]),int(cols[b]);use=(saved['idx']>=first)&(saved['idx']<=last)
                w=1/saved['var'][0,use];mean=float(np.sum(w*saved['d'][0,use])/w.sum())
                probability=fit['posterior'][a:b+1,s]
                # q labels component 1; the existing contrast is group 0 minus 1.
                segments.append(dict(first_bin=first,last_bin=last,n_test_bins=b-a+1,direction='lower' if s==2 else 'higher',
                    mean_diff=mean,z=float(mean*np.sqrt(w.sum())),posterior=float(probability.mean()),
                    min_bin_posterior=float(probability.min()),posterior_rule='every_retained_bin'))
            call.update(exact_segments=segments,log_evidence=fit['log_evidence'],states_enumerated=fit['states_enumerated'],
                        posterior=fit['posterior'].tolist(),reason='exact_conditional')
            if context=='depth_matched' and stage=='exact_initial':return old
            return [dict(old[0],segments=segments)]
        M.event_mixtures=events;M._score_candidate=scored;M.differences=difference
        try:yield
        finally:M.event_mixtures=event;M._score_candidate=score;M.differences=differences
