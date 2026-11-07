"""Read-only stage trace of membership, profile updates and merging on frozen inputs."""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import traceback

import numpy as np
import pandas as pd
from haplodose import clones, pipeline
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--assertion-source', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    arrays = args.out / 'arrays'
    arrays.mkdir()
    data = args.input_root / 'data'
    model_cells = args.input_root / 'model_cells.tsv'
    meta = pd.read_csv(model_cells, sep='\t')
    if any(c.startswith('truth_') for c in meta):
        raise ValueError('Truth columns present in inference metadata')
    s1 = set(np.flatnonzero(meta.specimen.eq('S1')))
    params = clones.CloneParams(split_mode='thin')
    record = start_run(args.out, dict(runner=Path(__file__), counts=data/'counts.h5',
        cells=model_cells, positions=data/'positions.tsv', assertions=args.assertion_source),
        asdict(params), {'scope':'Fixed development input, observational tracing only'}, mode='membership_trace')
    targets = {getattr(clones.SpecimenModel, x).__code__: x for x in
               ['shared_profiles', 'em', 'merge_unsupported', '_apply_splits', 'matched_differences']}
    targets[clones._fit_event_mixture.__code__] = 'event_fit'
    events, errors = [], []

    def observe(frame, event, result):
        if event != 'return' or frame.f_code not in targets:
            return
        try:
            loc = frame.f_locals
            cursor, m = frame, None
            while cursor is not None:
                if isinstance(cursor.f_locals.get('self'), clones.SpecimenModel):
                    m = cursor.f_locals['self']
                    break
                cursor = cursor.f_back
            if m is None or not set(m.rows.rows).issubset(s1):
                return
            stage = targets[frame.f_code]
            entry = dict(stage=stage, index=len(events), log_length=len(m.log))
            saved = {'rows': m.rows.rows.copy()}
            if stage == 'event_fit':
                outer = frame.f_back.f_locals
                if 'part' not in outer or 'cols' not in outer or result is None:
                    return
                cols = outer['cols'][outer['part']]
                # The fixed failing region, selected in the previous investigation.
                if cols[0] != 209 or cols[-1] != 217:
                    return
                saved.update({k: np.asarray(loc[k]).copy() for k in ['O','E','r']})
                saved.update(q=result['responsibility'].copy(), parameters=result['parameters'].copy(),
                    cells=outer['cells'].copy(), cols=cols.copy(),
                    Ytest=m.rows.take(outer['cells'], m.test_fold)[:,cols])
                entry['gain'] = float(result['gain'])
            else:
                for k in ['Yk','Lk','XLk','R','LL','pi','doses','labels']:
                    if k in loc and loc[k] is not None:
                        saved[k] = np.asarray(loc[k]).copy()
                if stage == 'shared_profiles':
                    saved['new_doses'] = result['dose'].copy()
                    saved['own_doses'] = result['dose_own'].copy() if 'dose_own' in result else result['dose'].copy()
                    entry['differences'] = loc.get('res', [])
                elif stage == 'em':
                    saved.update(R=result[2].copy(), LL=result[3].copy(), doses=np.asarray(result[0]), pi=result[1].copy())
                elif stage in ('merge_unsupported','_apply_splits'):
                    saved['new_doses'] = np.asarray(result[0]).copy()
                    entry['changed'] = bool(result[2])
                elif stage == 'matched_differences':
                    entry.update(pairs=loc['pairs'], differences=result)
            path = arrays / f'{len(events):04d}_{stage}.npz'
            np.savez_compressed(path, **saved)
            entry['arrays'] = str(path.relative_to(args.out))
            events.append(entry)
        except Exception:
            errors.append(traceback.format_exc())

    try:
        sys.setprofile(observe)
        try:
            pipeline.run(data/'counts.h5', model_cells, data/'positions.tsv', args.out/'fit',
                         params=params,chunk=97,counts_kind='raw',genome_build='hg38',
                         coordinate_system='one-based-closed',gene_namespace='symbol')
        finally:
            sys.setprofile(None)
        # Truth is joined only after the unmodified inference has finished.
        truth = pd.read_csv(data/'cells.tsv',sep='\t')
        tc = ['cell_id'] + [c for c in truth if c.startswith('truth_')]
        cells = pd.read_csv(args.out/'fit/cells.tsv.gz',sep='\t').merge(truth[tc],on='cell_id',validate='one_to_one')
        view=args.out/'assertion_view';view.mkdir()
        cells.to_csv(view/'cells.tsv.gz',sep='\t',index=False)
        for f in ['qc.json','differences.tsv']:shutil.copyfile(args.out/'fit'/f,view/f)
        sp=importlib.util.spec_from_file_location('original_assertions',args.assertion_source)
        mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
        gate='passed'
        try:mod.test_thin_records_mode_and_recovers_clones({'thin':view})
        except AssertionError:gate='failed'
        report=[]
        for e in events:
            z=np.load(args.out/e['arrays'])
            r=dict(stage=e['stage'],index=e['index'],log_length=e['log_length'])
            for k in ['doses','new_doses','own_doses']:
                if k in z:
                    r[k+'_region_means']=z[k][:,209:218].mean(axis=1).tolist()
                    r[k+'_region_spans']=np.ptp(z[k][:,209:218],axis=1).tolist()
            if 'R' in z:
                R=z['R'];r['effective_sizes']=R.sum(axis=0).tolist()
                r['mean_max_posterior']=float(R.max(axis=1).mean())
            if e['stage']=='event_fit':
                q=z['q'];g=q>.5;L=z['E'].sum(axis=1)
                # Proposed hard-pool purities under the fitted latent model.
                purity=[float(np.average(q[g==j],weights=L[g==j])) for j in [False,True]]
                r.update(gain=e['gain'],parameters=z['parameters'].tolist(),hard_pool_purities=purity,
                         attenuation=purity[1]-purity[0],q_quantiles=np.quantile(q,[0,.1,.25,.5,.75,.9,1]).tolist())
                actual=truth.iloc[z['rows'][z['cells']]].truth_clone.to_numpy()
                r['postfit_truth_table']=pd.crosstab(actual,g).to_dict()
            if 'differences' in e:
                r['difference_segments']=[d['segments'] for d in e['differences']]
            report.append(r)
        def default(x):
            if isinstance(x,np.generic):return x.item()
            raise TypeError(type(x).__name__)
        (args.out/'trace.json').write_text(json.dumps(events,indent=2,default=default)+'\n')
        (args.out/'summary.json').write_text(json.dumps(dict(gate=gate,events=report,errors=errors),indent=2,default=default)+'\n')
        if errors:raise RuntimeError('Trace errors; inspect summary')
        finish_run(args.out,record,backend='cpu')
    except Exception as error:
        finish_run(args.out,record,backend='cpu',error=error);raise


if __name__=='__main__':main()
