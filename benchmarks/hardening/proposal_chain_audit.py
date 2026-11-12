"""Independent post-fit identity, selection and NB-reference audit of proposals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.special import expit, logsumexp
from scipy.stats import nbinom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--trace-root',default='trace_v1')
    a=parser.parse_args();a.out.mkdir(exist_ok=False)
    run=a.task/a.trace_root
    assert json.loads((run/'report.json').read_text())['all_public_outputs_identical']
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',trace=run/'report.json'),{},
                  dict(scope='Truth used after fitting only;diagnostic purity not clone recovery'),mode='proposal_chain_audit')
    details=[];summaries=[];total_numeric=0;total_selections=0
    try:
        subprocess.run([sys.executable,'-I',str(Path(__file__).with_name('offset_recovery_audit.py')),
            '--out',str(a.out/'original_criteria'),'--task',str(a.task),'--arms',
            f'{a.trace_root}/bin/panel',f'{a.trace_root}/gene/panel'],check=True)
        for arm in ('bin','gene'):
            previous=a.task.parent/f'split_draw_stability_20261109/gate_v2/{arm}/default/panel'
            for seed in (20261031,20261107,20261108,20261109,20261110):
                left,right=[p/f'seed_{seed}/fit' for p in (run/arm/'panel',previous)]
                for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
                    pd.testing.assert_frame_equal(pd.read_csv(left/name,sep='\t'),pd.read_csv(right/name,sep='\t'),check_exact=True)
                assert np.array_equal(np.load(left/'posteriors.npy'),np.load(right/'posteriors.npy'),equal_nan=True)
            records=json.loads((run/arm/'trace.json').read_text());models={}
            for item in records:
                if item['kind']!='model':continue
                model=item['model'];seed=item['seed']
                z=np.load(run/arm/f'arrays/model_{model:02d}.npz')
                source=a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}'
                truth=pd.read_csv(source/'data/cells.tsv',sep='\t')
                metadata=pd.read_csv(source/'model_cells.tsv',sep='\t')
                assert truth.cell_id.tolist()==metadata.cell_id.tolist()
                rows=z['rows'];assert len(np.unique(rows))==len(rows)
                selected=truth.iloc[rows].reset_index(drop=True)
                assert selected.specimen.nunique()==1
                models[model]=dict(rows=rows,truth=selected,seed=seed,specimen=str(selected.specimen.iloc[0]),
                                   bins=pd.read_csv(item['bins'],sep='\t'))
            calls={};numeric=[]
            for item in records:
                if item['kind']=='event_fit':calls.setdefault(item['call'],[]).append(item)
            for call,events in calls.items():
                candidates=[r for r in records if r['kind']=='candidate' and r['call']==call and r['source'].startswith('event_')]
                eligible=[]
                for item in events:
                    z=np.load(item['arrays']);g=z['q']>.5
                    assert g.ndim==1 and len(g)==len(z['cells'])
                    if item['gain']>item['budget'] and 30<=int(g.sum())<=len(g)-30:
                        eligible.append((item,g))
                    if 'O' in z and 'parameters' in z:
                        O,E,r=z['O'],z['E'],z['r'];x=z['parameters']
                        ll=[nbinom.logpmf(O,r,r/(r+E*np.exp(eta))).sum(axis=1) for eta in x[:2]]
                        one=nbinom.logpmf(O,r,r/(r+E*np.exp(float(z['null_logtheta'])))).sum(axis=1)
                        joint=np.column_stack([ll[0]-np.logaddexp(0.,x[2]),ll[1]-np.logaddexp(0.,-x[2])])
                        mix=logsumexp(joint,axis=1);gain=float(np.sum(mix-one))
                        # Absolute SciPy log-PMF subtraction retains large common constants.
                        # Account for their cancellation; no scientific tolerance is changed.
                        budget=128*np.finfo(float).eps*(1+sum(float(np.abs(v).sum()) for v in (*ll,one,mix)))
                        assert abs(gain-item['gain'])<=budget+item['budget']
                        q=np.exp(joint[:,1]-mix)
                        assert np.max(np.abs(q-z['q']))<=budget
                        numeric.append(dict(model=item['model'],source=item['source'],gain=gain,
                                            difference=abs(gain-item['gain']),budget=budget))
                # Stable gain ordering and scalar sets independently reproduce proposal deduplication.
                eligible.sort(key=lambda pair:-pair[0]['gain']);kept=[];retained=[]
                for item,g in eligible:
                    minority=set(np.flatnonzero(~g if sum(g)>len(g)/2 else g).tolist())
                    if all(len(minority & old)/max(len(minority | old),1)<.7 for old in kept):
                        kept.append(minority);retained.append(item)
                    if len(retained)==5:break
                assert [r['source'] for r in candidates]==[r['source'] for r in retained]
                for returned,item in zip(candidates,retained):
                    candidate=np.load(returned['arrays']);window=np.load(item['arrays'])
                    assert np.array_equal(candidate['g'],window['q']>.5)
                    assert np.array_equal(candidate['cells'],window['cells'])
                total_selections+=1
            total_numeric+=len(numeric)
            write_json(numeric,a.out/f'{arm}_numerical.json')
            for item in records:
                if item['kind'] in ('model','candidate_call'):continue
                model=models[item['model']];z=np.load(item['arrays'])
                if item['kind']=='merge_assignment':
                    cells=np.arange(len(model['rows']));g=z['labels']
                    assert np.array_equal(g,z['assignment'].argmax(axis=1))
                else:cells=z['cells'];g=z['g'] if 'g' in z else (z['q']>.5).astype(int)
                truth=model['truth'].iloc[cells];t=truth.truth_clone.to_numpy()
                tumour=np.isin(t,['clone1','clone2']);represented=int(tumour.sum())
                table={}
                for label in sorted(set(g[tumour].tolist())):
                    table[str(label)]=[int(sum((g==label)&(t==name))) for name in ('clone1','clone2')]
                # Report the explicit contingency table, coverage and descriptive majority accuracy.
                correct=sum(max(pair) for pair in table.values())
                row=dict(arm=arm,seed=model['seed'],specimen=model['specimen'],kind=item['kind'],
                         source=item.get('source'),context=item.get('context'),model=item['model'],
                         arrays=item['arrays'],tumour_cells=represented,
                         total_true_tumour=int(model['truth'].truth_clone.isin(['clone1','clone2']).sum()),
                         contingency=table,majority_correct=correct,
                         majority_accuracy=correct/represented if represented else None,
                         accepted=item.get('accepted'),score=item.get('score'),members=item.get('members'),
                         call=item.get('call'),gain=item.get('gain'),budget=item.get('budget'))
                if 'cols' in z:
                    bins=model['bins'].iloc[z['cols']]
                    row['chromosomes']=sorted(set(bins.chromosome));row['start']=int(bins.start.min());row['end']=int(bins.end.max())
                if item['kind']=='agglomerate':
                    mapping={int(old):new for new,group in enumerate(item['members']) for old in group}
                    after=np.array([mapping[int(label)] for label in g])
                    final={str(label):[int(sum((after==label)&(t==name))) for name in ('clone1','clone2')] for label in sorted(set(after[tumour].tolist()))}
                    row.update(after_contingency=final,after_majority_accuracy=sum(max(v) for v in final.values())/represented if represented else None)
                details.append(row)
            for model_id,model in models.items():
                if model['specimen']!='S1':continue
                rows=[r for r in details if r['arm']==arm and r['model']==model_id]
                events=[r for r in rows if r['kind']=='event_fit' and r['chromosomes']==['chr12']]
                candidates=[r for r in rows if r['kind']=='candidate']
                scores=[r for r in rows if r['kind']=='score']
                summaries.append(dict(arm=arm,seed=model['seed'],
                    chr12_windows=[{k:v for k,v in r.items() if k not in ('arrays','contingency')} for r in events],
                    candidates=[{k:v for k,v in r.items() if k!='arrays'} for r in candidates],
                    confirmation=[{k:v for k,v in r.items() if k!='arrays'} for r in scores],
                    agglomeration=[r for r in rows if r['kind']=='agglomerate']))
        write_json(details,a.out/'details.json');write_json(summaries,a.out/'S1_chain.json')
        write_json(dict(all_checks_passed=True,numerical_chr12_windows=total_numeric,
                        independently_replayed_selections=total_selections,original_workflows=10,original_criteria_checked=60,
                        production_changed=False,biological_specimens=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU independent SciPy NB, scalar-set selection and post-fit truth joins')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
