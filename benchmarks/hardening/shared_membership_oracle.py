"""Postfit diagnostic of membership-information loss; oracle is never deployable."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.clones import CloneParams,_thin_diff_params
from haplodose.profile import DIFF_LEVELS,DIFF_PRIOR,credible_difference_runs
from haplodose.provenance import start_run,finish_run


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--trace-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);path=Path(__file__).with_name('shared_membership_hmm.py')
    spec=importlib.util.spec_from_file_location('exact_oracle_model',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    inputs={'runner':__file__,'model':path}
    for seed in (20261108,20261109,20261110):
        inputs[f'{seed}:truth']=a.panel_root/f'seed_{seed}/data/cells.tsv'
        directory=a.trace_root/f'trace_{seed}'
        for f in ('trace.json','fit/bins.tsv','fit/specimens.tsv'):inputs[f'{seed}:{f}']=directory/f
        for event in json.loads((directory/'trace.json').read_text()):
            if event['stage']=='event_fit':inputs[f'{seed}:event{event["index"]}']=directory/event['arrays']
    rec=start_run(a.out,inputs,{},dict(role='postfit membership attribution',oracle_is_deployable=False,
        expectation='If oracle membership alone repairs support, classification information is limiting; otherwise exposure/state mismatch remains'),mode='exact_membership_oracle')
    try:
        rows=[];pp=_thin_diff_params(CloneParams(split_mode='thin').profile)
        for seed in (20261108,20261109,20261110):
            meta=pd.read_csv(inputs[f'{seed}:truth'],sep='\t');directory=a.trace_root/f'trace_{seed}'
            bins=pd.read_csv(directory/'fit/bins.tsv',sep='\t');chrom=bins.chromosome.to_numpy()
            beta=float(pd.read_csv(directory/'fit/specimens.tsv',sep='\t').set_index('specimen').loc['S1','beta'])
            for event in json.loads((directory/'trace.json').read_text()):
                if event['stage']!='event_fit':continue
                z=np.load(directory/event['arrays']);cols=z['cols'];global_rows=z['rows'][z['cells']]
                selected=meta.iloc[global_rows]
                if not selected.specimen.eq('S1').all():raise ValueError('Specimen mismatch')
                qtruth=selected.truth_clone.eq('clone2').to_numpy(float)
                for label,q in [('soft',z['q']),('hard',(z['q']>.5).astype(float)),('oracle_truth',qtruth)]:
                    fitted=m.build(z['Ytest'],z['E'],z['r'],q,DIFF_LEVELS*beta*np.log(1.5),DIFF_PRIOR,pp.p_change,
                        left_neutral=bool(cols[0]>0 and chrom[cols[0]-1]==chrom[cols[0]]),
                        right_neutral=bool(cols[-1]+1<len(chrom) and chrom[cols[-1]+1]==chrom[cols[-1]]))
                    fit=m.enumerate_states(fitted);reset=np.r_[True,np.zeros(len(cols)-1,bool)]
                    cores=[(int(cols[x]),int(cols[y]),int(s)) for x,y,s in credible_difference_runs(fit['path'],fit['posterior'],reset,pp)]
                    row=dict(seed=seed,membership=label,log_evidence=fit['log_evidence'],cores=cores,
                        posterior=fit['posterior'].tolist(),beta=beta,cell_ids=selected.cell_id.tolist(),
                        interpretation='Fitted exposure/dispersion and original state grid retained; mean centering follows supplied memberships; not recoverable oracle performance')
                    rows.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('cell_ids','posterior')}),flush=True)
        (a.out/'summary.json').write_text(json.dumps(rows,indent=2));finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
