"""Join audited joint draws and report the fixed 80-family with numerical bounds."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta,binom
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    import importlib.util
    spec=importlib.util.spec_from_file_location('uncertainty_summary_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def verified(directory):
    manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--chunks',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    candidate=load('bulk_amplitude_uncertainty');stats=load('statistics');reference=load('gene_factorial_precision_audit')
    chunks=json.loads(a.chunks.read_text());assert len(chunks)==20
    inputs=dict(script=Path(__file__),candidate=Path(candidate.__file__),statistics=Path(stats.__file__),
                independent_arithmetic=Path(reference.__file__),chunks=a.chunks,design=a.task/'DESIGN.md')
    cursor=0;origin=None;definition=None
    for i,chunk in enumerate(chunks):
        directory=Path(chunk['path']);audit=Path(chunk['audit']);m=verified(directory);verified(audit)
        check=json.loads((audit/'report.json').read_text());assert check['all_checks_passed']
        report=json.loads((directory/'report.json').read_text())
        assert report['start']==cursor==chunk['start'] and report['stop']==chunk['stop'];cursor=chunk['stop']
        current=json.loads((directory/'definition.json').read_text())
        if definition is None:definition=current;origin=m['input']['files']
        else:assert current==definition
        for key in ('script','design','profile',*[f'windows_{j}' for j in range(9)],*[f'scores_{j}' for j in range(9)]):
            assert m['input']['files'][key]['sha256']==origin[key]['sha256']
        for name in ('run_manifest.json','draws.npy','paired.npy','indices.npy','thresholds.npy','definition.json'):
            inputs[f'chunk_{i}:{name}']=directory/name
        inputs[f'audit_{i}']=audit/'report.json'
    assert cursor==candidate.TOTAL and file_record(inputs['candidate'])['sha256']==origin['script']['sha256']
    bank,corr,methods,bank_inputs=candidate.read_bank(a.task.parents[1])
    inputs.update({'bank_'+key:value for key,value in bank_inputs.items()})
    rec=start_run(a.out,inputs,dict(total=cursor,scientific_family=80,monte_carlo_family=240),
        dict(scope='Approximate sampling uncertainty conditional on exposed patients and finite reference bank; no validated biological coverage'),mode='bulk_amplitude_uncertainty_summary')
    try:
        draws=np.lib.format.open_memmap(a.out/'draws.npy',mode='w+',dtype=float,shape=(cursor,80))
        paired=np.lib.format.open_memmap(a.out/'paired.npy',mode='w+',dtype=np.uint8,shape=(cursor,80))
        for i,chunk in enumerate(chunks):
            lo,hi=chunk['start'],chunk['stop']
            values=np.load(inputs[f'chunk_{i}:draws.npy'],mmap_mode='r');counts=np.load(inputs[f'chunk_{i}:paired.npy'],mmap_mode='r')
            assert values.shape==counts.shape==(hi-lo,80)
            assert np.array_equal(np.isfinite(values),counts>0)
            draws[lo:hi]=values;paired[lo:hi]=counts
        draws.flush();paired.flush();rows=[];mc_tail=.05/(2*240)
        w=bank[0];ci=np.flatnonzero(w['split']=='calibration');ei=np.flatnonzero(w['split']=='evaluation')
        point={}
        for m in methods:
            threshold=candidate.Calibration(w[m],w['neutral']&w['common'],ci).threshold(ci)
            score=candidate.LossCurve(w[m],w['loss'],w['neutral'],w['gain'],w['common']).at(threshold)
            point[m]=np.column_stack([corr[0][m],score])
        for j,contrast in enumerate(definition['contrasts']):
            delta=definition['point_effects'][j];values=np.asarray(draws[:,j]);finite=np.isfinite(values)
            b=values[finite];n=len(b);undefined=cursor-n
            endpoint=candidate.ENDPOINTS.index(contrast['endpoint'])
            differences=(point[contrast['method']][:,endpoint]-point[contrast['control']][:,endpoint])[ei]
            differences=differences[np.isfinite(differences)]
            assert len(differences)==definition['point_paired_patients'][j]
            estimable=bool(np.isfinite(delta) and len(differences)>=2 and np.ptp(differences)>0 and undefined==0 and np.ptp(b)>0)
            row=dict(**contrast,mean_diff_oriented=delta,mean_diff_raw=delta*contrast['direction'],
                point_paired_patients=definition['point_paired_patients'][j],valid_draws=n,undefined_draws=undefined,
                undefined_fraction=undefined/cursor,undefined_probability_upper_marginal95=float(beta.isf(.05,undefined+1,cursor-undefined)) if undefined<cursor else 1.,
                minimum_paired_patients=int(paired[:,j].min()),maximum_paired_patients=int(paired[:,j].max()),
                testable=estimable,p=1.,p_mc_lower=0.,p_mc_upper=1.,family_ci_low=np.nan,family_ci_high=np.nan,
                family_ci_low_mc=[-np.inf,np.inf],family_ci_high_mc=[-np.inf,np.inf])
            if n and np.isfinite(delta):
                errors=b-delta;ordered=np.sort(errors)
                # Conditional defined-draw summaries remain even when no unconditional gate is available.
                for key,q in (('ci_low',.975),('ci_high',.025),('family_ci_low',1-.05/160),('family_ci_high',.05/160)):
                    value=float(delta-np.quantile(errors,q));independent,guard=reference.basic_endpoint(ordered,delta,q)
                    assert abs(value-independent)<=guard,(j,key,value,independent,guard)
                    row[key]=value
                if estimable:
                    k=int(np.count_nonzero(np.abs(errors)>=abs(delta)));row['exceedances']=k
                    row['p']=float(Fraction(k+1,n+1))
                    row['p_mc_lower']=0. if k==0 else float(beta.ppf(mc_tail,k,n-k+1))
                    row['p_mc_upper']=1. if k==n else float(beta.isf(mc_tail,k+1,n-k))
                    for key,q in (('family_ci_low',1-.05/160),('family_ci_high',.05/160)):
                        lo=int(binom.ppf(mc_tail,n,q));hi=int(binom.isf(mc_tail,n,q))+1
                        assert binom.cdf(lo-1,n,q)<=mc_tail and binom.sf(hi-1,n,q)<=mc_tail
                        lower=-np.inf if lo==0 else float(ordered[lo-1]);upper=np.inf if hi>n else float(ordered[hi-1])
                        row[key+'_mc']=[delta-upper,delta-lower]
            rows.append(row)
        # The entire 80-family remains, including untestable entries fixed at p=1.
        for key in ('p','p_mc_lower','p_mc_upper'):
            adjusted=stats.holm([r[key] if r['testable'] else 1. for r in rows])
            independent=reference.holm([r[key] if r['testable'] else 1. for r in rows])
            assert np.array_equal(adjusted,np.asarray(independent))
            for row,value in zip(rows,adjusted):row[key+'_holm']=float(value)
        for row in rows:
            row['conditional_gate_passed']=bool(row['testable'] and row['p_mc_upper_holm']<.05 and row['family_ci_low_mc'][0]>0)
            row['numerically_unresolved']=bool(row['testable'] and (
                row['p_mc_lower_holm']<=.05<=row['p_mc_upper_holm'] or row['family_ci_low_mc'][0]<=0<=row['family_ci_low_mc'][1]))
        summary=dict(draws=cursor,scientific_family=80,monte_carlo_family=240,
            comparisons=rows,conditional_gates_passed=sum(r['conditional_gate_passed'] for r in rows),
            numerically_unresolved=sum(r['numerically_unresolved'] for r in rows),
            all_chunk_audits_passed=True,normal_reference_bank_draws=8,comparison_patients=24,calibration_patients=24,
            interval_scope='Approximate joint bootstrap conditional on eight frozen normal-reference resamples and exposed patients; empirical coverage unvalidated',
            precision_scope='Paired callers only,not unconditional cohort precision',method_promoted=False,biological_coverage_validated=False)
        (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(a.out,rec,backend='CPU rational interpolation,independent Holm and exact binomial Monte Carlo bounds')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
