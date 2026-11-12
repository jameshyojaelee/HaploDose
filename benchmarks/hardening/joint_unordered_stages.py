"""Bounded joint-count recovery diagnostic, never independent confirmation.

Returned short event candidates retain the discovery search/ranking but use an
unanchored full-count mixture to assign cells and assess unordered trajectories.
Long/non-event candidates retain the old procedure, explicitly logged. Selection
and nuisance fitting are adaptive; no posterior calibration claim is made.
"""
from contextlib import contextmanager
import importlib.util
from pathlib import Path
import sys
import numpy as np
from scipy.special import expit
from haplodose import clones, profile

def load(name):
    spec=importlib.util.spec_from_file_location('joint_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


balanced_subset=load('shared_membership_stages').balanced_subset


@contextmanager
def joint_stages(merge, trace):
    if merge not in {'hard', 'margin'}:raise ValueError(merge)
    model=load('shared_membership_hmm');unordered=load('label_invariant_membership')
    with load('membership_weighted_stages').weighted_stages('merge' if merge=='margin' else 'hard',trace,
            record_memberships=True):
        if sys.getprofile() is not None:raise RuntimeError('Existing profiler')
        M=clones.SpecimenModel;event,score,differences=M.event_mixtures,M._score_candidate,M.differences
        cache={};active=None

        def fit(self,item,cells,choose):
            cols=item['cols'];p=item['p'];share=self.fold_scale['prop']
            if not np.all(np.diff(cols)==1) or np.any(self.reset[cols[1:]]):raise ValueError('Discontinuous event')
            Y=self.rows.take(cells)[:,cols]
            built=model.build(Y[choose],item['E'][choose]/share,item['r'][choose]/share,
                np.full(len(choose),p),profile.DIFF_LEVELS*self.beta*np.log(1.5),profile.DIFF_PRIOR,
                self.diff_pp.p_change,left_neutral=bool(not self.reset[cols[0]]),
                right_neutral=bool(cols[-1]+1<self.J and not self.reset[cols[-1]+1]))
            result=unordered.enumerate_unordered(built)
            nops=64+4*Y[choose].size;eps=np.finfo(float).eps
            allowance=4*nops*eps/(1-nops*eps)*max(1.,float(np.abs(built['ll']).sum()))
            support=result['posterior'][np.arange(len(cols)),result['path']]
            if np.min(np.abs(support-self.diff_pp.posterior))<=allowance:
                raise ArithmeticError('State decision overlaps numerical allowance')
            result['numerical_allowance']=allowance
            return result

        def events(self,*args,**kwargs):
            if self.alleles is not None or self.graph is not None or self.params.split_mode!='thin' or self.params.split_null_draws:
                raise ValueError('Expression-only nonspatial research procedure without selected-parent null')
            cells=np.asarray(kwargs.get('cells',args[5] if len(args)>5 else None))
            if cells.ndim!=1:raise ValueError('Missing cell identities')
            captured={}
            def observe(frame,kind,value):
                if kind!='return' or frame.f_code is not clones._fit_event_mixture.__code__ or not isinstance(value,dict):return
                loc=frame.f_back.f_locals
                if loc.get('self') is not self or 'cols' not in loc or 'part' not in loc:return
                cols=loc['cols'][loc['part']]
                source=f"event_chrom{loc['c']}_bins{cols[0]}-{cols[-1]}"
                captured[source]=dict(q=value['responsibility'].copy(),p=float(expit(value['parameters'][2])),
                    E=loc['E'].copy(),r=loc['r_nb'].copy(),cols=cols.copy())
            sys.setprofile(observe)
            try:found=event(self,*args,**kwargs)
            finally:sys.setprofile(None)
            key=(id(self),cells.tobytes());cache[key]={};returned=[]
            for g,source in found:
                item=captured[source]
                if not np.array_equal(item['q']>.5,g.astype(bool)):raise ValueError('Original probability alignment')
                entry=dict(stage='joint_proposal',source=source,specimen_rows=self.rows.rows.tolist(),
                    proposal_cells=cells.tolist(),original_labels=g.tolist(),confirmation_status='joint_not_independent')
                trace.append(entry)
                if 3**len(item['cols'])>1000000:
                    entry['fallback']='declared_exact_state_cap';returned.append((g,source));continue
                fitted=fit(self,item,cells,np.arange(len(cells)))
                labels=(fitted['membership']>.5).astype(int)
                entry.update(membership=fitted['membership'].tolist(),labels=labels.tolist(),
                    neutral_probability=fitted['neutral_trajectory_probability'],p=item['p'],
                    numerical_allowance=fitted['numerical_allowance'],posterior=fitted['posterior'].tolist(),
                    path=fitted['path'].tolist())
                cores=profile.credible_difference_runs(fitted['path'],fitted['posterior'],
                    np.r_[True,np.zeros(len(item['cols'])-1,bool)],self.diff_pp)
                if not cores:
                    entry['rejected']='unchanged_joint_core_rule';continue
                if min(np.bincount(labels,minlength=2))<self.params.min_cells:
                    entry['rejected']='unchanged_minimum_cells';continue
                if np.min(np.abs(fitted['membership']-.5))<=fitted['numerical_allowance']:
                    entry['rejected']='unresolved_membership_numerics';continue
                item['fit']=fitted;cache[key][source]=item;returned.append((labels,source))
            return returned

        def scored(self,k,cells,g,source):
            nonlocal active
            item=cache.get((id(self),np.asarray(cells).tobytes()),{}).get(source)
            if item is None:return score(self,k,cells,g,source)
            entry=dict(stage='joint_scoring',source=source,parent=int(k),calls=[],
                specimen_rows=self.rows.rows.tolist(),proposal_cells=np.asarray(cells).tolist())
            trace.append(entry);before=active;active=(self,item,np.asarray(cells),g,entry)
            try:result=score(self,k,cells,g,source)
            finally:active=before
            entry['accepted']=result is not None
            return result

        def difference(self,*args,**kwargs):
            if active is None or active[0] is not self:return differences(self,*args,**kwargs)
            _,item,cells,g,entry=active;caller=sys._getframe(1)
            if caller.f_code is M.agglomerate.__code__:
                if caller.f_locals['A']!=[0] or caller.f_locals['B']!=[1]:raise ValueError('Unexpected binary pair')
                context='initial'
            elif caller.f_code.co_name=='matched_differences' and 'W' in caller.f_locals:context='depth_matched'
            else:raise ValueError('Unknown scoring context')
            saved={}
            def observe(frame,kind,value):
                if kind=='return' and frame.f_code is profile.difference_segments.__code__:
                    saved.update({k:frame.f_locals[k].copy() for k in ('d','var','idx')})
            sys.setprofile(observe)
            try:old=differences(self,*args,**kwargs)
            finally:sys.setprofile(None)
            if len(old)!=1 or saved['d'].shape[0]!=1:raise ValueError('Expected one contrast')
            choose=np.arange(len(cells));sizes=None
            if context=='depth_matched':choose,sizes=balanced_subset(self.x[cells],g,self.rows.rows[cells])
            call=dict(context=context,original_segments=old[0]['segments'],selected_cells=choose.tolist(),
                per_stratum_per_group=sizes,selection_status='joint-data labels; selection not conditioned out')
            entry['calls'].append(call)
            if len(choose)<2*self.params.min_cells:
                call.update(joint_segments=[],reason='insufficient_depth_overlap')
                return [dict(old[0],segments=[],testable=False,untestable_reason='insufficient_depth_overlap')]
            fitted=item['fit'] if context=='initial' else fit(self,item,cells,choose)
            # Subset labels may use a different first nonneutral locus. Align
            # ONE global permutation to the proposal, never individual bins.
            overlap=float(np.sum((2*g[choose]-1)*(2*fitted['membership']-1)))
            if abs(overlap)<=len(choose)*fitted['numerical_allowance']:
                raise ArithmeticError('Component-to-group orientation unresolved')
            if overlap<0:
                fitted=dict(fitted,path=2-fitted['path'],posterior=fitted['posterior'][:,::-1])
            call['global_label_exchange_to_proposal']=bool(overlap<0)
            cols=item['cols'];segments=[]
            for a,b,s in profile.credible_difference_runs(fitted['path'],fitted['posterior'],
                    np.r_[True,np.zeros(len(cols)-1,bool)],self.diff_pp):
                first,last=int(cols[a]),int(cols[b]);use=(saved['idx']>=first)&(saved['idx']<=last)
                w=1/saved['var'][0,use];mean=float(np.sum(w*saved['d'][0,use])/w.sum())
                probability=fitted['posterior'][a:b+1,s]
                segments.append(dict(first_bin=first,last_bin=last,n_test_bins=b-a+1,
                    direction='lower' if s==2 else 'higher',mean_diff=mean,z=float(mean*np.sqrt(w.sum())),
                    posterior=float(probability.mean()),min_bin_posterior=float(probability.min()),
                    posterior_rule='every_retained_bin',estimand='unordered_joint_RNA_trajectories'))
            call.update(joint_segments=segments,log_evidence=fitted['log_evidence'],
                posterior=fitted['posterior'].tolist(),path=fitted['path'].tolist(),
                neutral_probability=fitted['neutral_trajectory_probability'],numerical_allowance=fitted['numerical_allowance'])
            return [dict(old[0],segments=segments)]

        M.event_mixtures=events;M._score_candidate=scored;M.differences=difference
        try:yield
        finally:M.event_mixtures=event;M._score_candidate=score;M.differences=differences
