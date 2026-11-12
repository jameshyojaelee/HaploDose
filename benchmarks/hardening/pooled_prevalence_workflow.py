"""Full original recovery with a discovery-locked pooled-prevalence safeguard."""
import argparse
from contextlib import contextmanager
from functools import lru_cache
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('pooled_workflow_'+name,
        Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


@contextmanager
def safeguard(method,registry,mode,trace,directory,records,guard_directory):
    """Fit during discovery preparation, before the original confirmation call."""
    original_load=method.load
    fitted={}
    def loader(name):
        module=original_load(name)
        if name!='depth_conditioned_genes':return module
        original_prepare=module.prepare
        def prepare(*args,**kwargs):
            parameters=original_prepare(*args,**kwargs)
            frame=sys._getframe(1)
            if frame.f_code.co_name!='scored' or frame.f_globals is not method.__dict__:
                raise RuntimeError('Discovery preparation called outside the fixed adapter')
            loc=frame.f_locals
            self,item,info,cells=loc['self'],loc['item'],loc['info'],loc['cells']
            source=loc['source']
            key=(id(self.rows),np.asarray(cells).tobytes(),source)
            number=len(records)
            folder=guard_directory/f'proposal_{number:05d}'
            folder.mkdir()
            record=dict(number=number,case_seed=info['case_seed'],source=source,
                folder=str(folder),status='discovery_fitting',confirmation_counts_opened_during_guard_fit=False)
            records.append(record)
            write_json(records,guard_directory/'records.json')
            if not parameters['identified']:
                record.update(status='unidentified_discovery_support')
                fitted[key]=dict(record=record,available=False)
                write_json(records,guard_directory/'records.json')
                return parameters
            choose=np.asarray(parameters['selected'],int)
            strata=np.asarray(parameters['strata'],int)[choose]
            offsets=np.exp(np.asarray(parameters['offsets'])[choose,None])
            share=info['share']
            # Construct test-fold EXPECTATIONS from discovery inputs only,
            # in the exact offline-screen operation order. No C counts read.
            array_module=load('gene_event_proposal')
            expected_bins=item['E']*(self.fold_scale['test']/self.fold_scale['prop'])
            _,Eg,rg,genes,local=array_module.region_arrays(loc['X'],info['gbin'],
                loc['fraction'],info['resolved_gene_size']*(1-share),item['cols'],expected_bins)
            np.testing.assert_array_equal(genes,loc['genes'])
            D=np.asarray(args[0])[choose]
            EC=Eg[choose]*offsets
            rC=np.broadcast_to(rg,Eg.shape)[choose]
            ED,rD=EC*share/(1-share),rC*share/(1-share)
            q=np.asarray(item['q'])[choose]
            np.savez_compressed(folder/'discovery.npz',D=D,ED=ED,rD=rD,EC=EC,rC=rC,
                strata=strata,q=q,selected=choose,genes=genes,global_rows=self.rows.rows[cells][choose])
            models=[]
            try:
                for order in (24,48,96):
                    print(json.dumps(dict(case=info['case_seed'],proposal=number,order=order,status='guard_discovery_fit')),flush=True)
                    old=load('depth_state_predictive').fit(D,ED,rD,strata,q,order)
                    models.append(load('profiled_prevalence').fit(D,ED,rD,strata,q,old,shared=True))
            except ArithmeticError as error:
                record.update(status='unresolved_discovery_optimizer',reason=str(error))
                fitted[key]=dict(record=record,available=False)
                write_json(dict(models=models,record=record,confirmation_counts_opened=False),folder/'DISCOVERY_MODELS_LOCKED.json')
                write_json(records,guard_directory/'records.json')
                return parameters
            write_json(dict(models=models,record=record,confirmation_counts_opened=False),folder/'DISCOVERY_MODELS_LOCKED.json')
            fitted[key]=dict(record=record,available=True,models=models,D=D,ED=ED,rD=rD,EC=EC,rC=rC,
                strata=strata,choose=choose,genes=genes,cells=np.asarray(cells).copy(),item=item,info=info)
            record.update(status='discovery_models_locked')
            write_json(records,guard_directory/'records.json')
            return parameters
        module.prepare=prepare
        return module
    method.load=loader
    try:
        with method.confirmation(registry,mode,trace,directory):
            original=clones.SpecimenModel._score_candidate
            def score(self,k,cells,g,source):
                result=original(self,k,cells,g,source)
                if result is None:return None
                key=(id(self.rows),np.asarray(cells).tobytes(),source)
                if key not in fitted:
                    # Retain the existing non-event/declared long-event fallback
                    # mechanism. Whole-workflow calibration must cover it too.
                    return result
                state=fitted[key];record=state['record']
                if not state['available']:
                    record['guard_passed']=False
                    write_json(records,guard_directory/'records.json')
                    return None
                observation=state['item']['conditioned_gene_observation']
                choose=state['choose']
                C=np.asarray(observation['counts'])[choose]
                offsets=np.exp(np.asarray(state['item']['conditioned_parameters']['offsets'])[choose,None])
                EC=np.asarray(observation['expected'])[choose]*offsets
                rC=np.broadcast_to(observation['size'],observation['counts'].shape)[choose]
                np.testing.assert_array_equal(EC,state['EC'])
                np.testing.assert_array_equal(rC,state['rC'])
                full=state['info']['X'][state['cells']][:,state['genes']].toarray()[choose]
                np.testing.assert_array_equal(full,state['D']+C)
                rows=[]
                folder=Path(record['folder'])
                for model in state['models']:
                    value=load('profiled_prevalence').predict(state['D'],C,state['ED'],EC,
                        state['rD'],rC,state['strata'],model)
                    arrays={name:value.pop(name)for name in('per_cell','discovery_membership','joint_membership')}
                    np.savez_compressed(folder/f"prediction_{model['order']}.npz",**arrays)
                    rows.append(dict(order=model['order'],**value))
                diagnostic=2*max(abs(rows[i]['score']-rows[i-1]['score'])for i in(1,2))+sum(r['rounding_budget']for r in rows)
                resolved=bool(abs(rows[-1]['score'])>diagnostic)
                passed=bool(resolved and rows[-1]['score']>0)
                record.update(status='confirmation_predicted',quadrature=rows,
                    score=rows[-1]['score'],numerical_uncertainty_diagnostic=diagnostic,
                    resolved=resolved,guard_passed=passed,fold_count_conservation_exact=True,
                    discovery_fit_reads_confirmation=False)
                write_json(record,folder/'report.json')
                write_json(records,guard_directory/'records.json')
                print(json.dumps(dict(case=record['case_seed'],proposal=record['number'],
                    status='guard_result',score=record['score'],passed=passed)),flush=True)
                return result if passed else None
            clones.SpecimenModel._score_candidate=score
            try:yield
            finally:clones.SpecimenModel._score_candidate=original
    finally:method.load=original_load


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--control',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    guards=args.out/'guards';guards.mkdir()
    records=[]
    prior=args.task.parent/'adjacent_reference_noise_20261110/candidate_v2'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',frozen_screen=args.task.parent/'pooled_prevalence_20261110/summary_v1/report.json',
        previous_recovery=prior/'report.json',method=Path(__file__).with_name('gene_shared_confirmation.py'),
        runner=Path(__file__).with_name('adjacent_reference_noise.py'),optimizer=Path(__file__).with_name('profiled_prevalence.py'),
        fixed_fit=Path(__file__).with_name('depth_state_predictive.py'),integrals=Path(__file__).with_name('state_predictive_confirmation.py'))
    inputs['runtime_control']=args.task/'runtime_control_v1/report.json'
    assert json.loads(inputs['runtime_control'].read_text())['all24arrays_exact']
    if not args.control:
        inputs['strict_workflow_control']=args.task/'strict_control_v1/report.json'
        assert json.loads(inputs['strict_workflow_control'].read_text())['strict_control']['exact']
    assert json.loads(inputs['frozen_screen'].read_text())['screen_passed']
    rec=start_run(args.out,inputs,dict(control=args.control,guard_mode='shared_prevalence_discovery_locked'),
        dict(scope='Original full recovery, including merge; not a specimen-error or biological claim'),mode='pooled_prevalence_workflow')
    adjacent=load('adjacent_reference_noise');original_load=adjacent.load;argv=sys.argv
    def loader(name):
        runner=original_load(name)
        if name!='gene_shared_recovery' or args.control:return runner
        old_load=runner.load
        def select(name):
            method=old_load(name)
            if name!='gene_shared_confirmation':return method
            original_confirmation=method.confirmation
            @contextmanager
            def context(registry,mode,trace,directory):
                # safeguard calls this original context while its preparation
                # loader is installed; avoid wrapping it recursively.
                method.confirmation=original_confirmation
                try:
                    with safeguard(method,registry,mode,trace,directory,records,guards):yield
                finally:method.confirmation=context
            method.confirmation=context
            return method
        runner.load=select
        return runner
    adjacent.load=loader
    try:
        sys.argv=[str(adjacent.__file__),'--task',str(args.task),'--out',str(args.out/'full')]
        adjacent.main()
        report=json.loads((args.out/'full/report.json').read_text())
        exact=None
        if args.control:
            before=json.loads(inputs['previous_recovery'].read_text())
            assert before['recovered']==report['recovered']==5
            assert [r['criteria']for r in before['cases']]==[r['criteria']for r in report['cases']]
            arrays,tables,folds=0,0,0
            for row in report['cases']:
                a=prior/'workflow/candidate/panel'/f"seed_{row['seed']}"/'fit'
                b=args.out/'full/workflow/candidate/panel'/f"seed_{row['seed']}"/'fit'
                for file in a.iterdir():
                    if file.suffix=='.npy':np.testing.assert_array_equal(np.load(file),np.load(b/file.name));arrays+=1
                    elif file.name.endswith(('.tsv','.tsv.gz')):
                        pd.testing.assert_frame_equal(pd.read_csv(file,sep='\t'),pd.read_csv(b/file.name,sep='\t'),check_exact=True);tables+=1
            for file in(prior/'workflow/candidate/folds').glob('*.npz'):
                with np.load(file)as a,np.load(args.out/'full/workflow/candidate/folds'/file.name)as b:
                    assert a.files==b.files
                    for key in a.files:np.testing.assert_array_equal(a[key],b[key])
                folds+=1
            assert folds==15
            exact=dict(criteria=30,arrays=arrays,tables=tables,folds=folds,exact=True)
        write_json(dict(control=args.control,recovered=report['recovered'],cases=report['cases'],
            strict_control=exact,guard_preparations=len(records),
            guard_predictions=sum(r['status']=='confirmation_predicted'for r in records),
            guard_unresolved=sum(r['status']=='unresolved_discovery_optimizer'for r in records),
            full_workflow_including_merge=True,false_split_control_established=False,
            production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        write_json(records,guards/'records.json')
        finish_run(args.out,rec,backend='CPU complete raw-count workflow with pre-confirmation shared-prevalence model fits')
    except BaseException as error:
        write_json(records,guards/'records.json');finish_run(args.out,rec,error=error);raise
    finally:adjacent.load,sys.argv=original_load,argv


if __name__=='__main__':main()
