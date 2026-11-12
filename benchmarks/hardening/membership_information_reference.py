"""Oracle membership information under the frozen state-plus-CNA generator.

An unavailable noise-free amplitude is an information reference, not a deployable
predictor, biological truth model, or per-realization accuracy bound. Cell IDs
are assumed to convey no genotype information; the generator's row order cannot
be used as a classifier. Fixed balanced class sizes are additionally supplied to
the top-half oracle, so it has more information than the inference procedure.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linear_sum_assignment
from scipy.stats import norm
from sklearn.metrics import adjusted_rand_score
from haplodose.provenance import start_run,finish_run,file_record


def accuracy(truth,predicted,eligible):
    """One-to-one matching; unassigned/reference-like/extra groups receive no credit."""
    classes=np.unique(truth);eligible=np.asarray(eligible)
    if not len(eligible):return 0.
    counts=np.array([[np.count_nonzero((truth==a)&(predicted==b)) for b in eligible] for a in classes])
    row,col=linear_sum_assignment(-counts)
    return float(counts[row,col].sum()/len(truth))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--rank-source',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    # Minimal endpoint counterexamples, including a split into four pure groups.
    truth=np.array([0,0,1,1])
    assert accuracy(truth,np.array([7,7,9,9]),[7,9])==1.
    assert accuracy(truth,np.array([1,2,3,4]),[1,2,3,4])==.5
    assert accuracy(truth,np.zeros(4,int),[])==0.
    assert accuracy(truth,np.array([7,7,0,0]),[7])==.5
    generator=a.rank_source/'global_rank_generator.py'
    design=json.loads((a.rank_source/'prespecified_design.json').read_text())
    if hashlib.sha256(generator.read_bytes()).hexdigest()!=design['generator_sha256']:
        raise ValueError('Frozen observation model changed')
    audit=a.task/'null_audit_v1/report.json'
    if not audit.is_file():raise ValueError('Full workflow audit is required first')
    inputs={'runner':__file__,'generator':generator,'design':a.rank_source/'prespecified_design.json',
            'audit':audit,'summary':a.task/'null_summary_v1/summary.json','evaluation':a.task/'null_summary_v1/evaluation.json'}
    for stage in ('merge_control','exact_matched'):
        inputs[stage]=a.task/'null_shards_v1'/stage/'eval_two/records.jsonl'
    delta=.6*np.log(1.5);sigma=.25
    analytic=float(norm.cdf(delta/(2*sigma)))
    a0,e0=quad(lambda x:.5*norm.pdf(x,loc=0,scale=sigma),-np.inf,delta/2)
    a1,e1=quad(lambda x:.5*norm.pdf(x,loc=delta,scale=sigma),delta/2,np.inf)
    if abs(analytic-a0-a1)>e0+e1+8*np.finfo(float).eps:raise ArithmeticError('Gaussian Bayes-risk reference disagreement')
    rec=start_run(a.out,inputs,dict(delta=delta,state_sd=sigma),
        dict(scope='Postfit conditional simulation information diagnostic; not biological validation',
             oracle='Noise-free total regional amplitude; true generator parameters; optional exact balanced class sizes',
             assumptions='Exchangeable cell identities; given amplitude and depth, count law carries no further genotype information'),mode='membership_information_reference')
    try:
        evaluations=json.loads(inputs['evaluation'].read_text());rows=[];fingerprints={}
        for stage in ('merge_control','exact_matched'):
            with inputs[stage].open() as stream:
                for i,line in enumerate(stream):
                    record=json.loads(line)
                    if record['seed_key']!=[20261107,20,2,i]:raise ValueError('Evaluation identity/order changed')
                    array_record=record['assignments']
                    if file_record(array_record['path'])['sha256']!=array_record['sha256']:raise ValueError('Assignment artifact changed')
                    with np.load(array_record['path']) as z:
                        truth=z['simulated_CNA_labels'];labels=z['raw_MAP'];latent=z['latent_state'];eligible=z['eligible_tumour_groups']
                    if not np.array_equal(np.bincount(truth),[300,300]):raise ValueError('Fixed alternative changed')
                    key=hashlib.sha256(truth.tobytes()+latent.tobytes()).hexdigest()
                    if stage=='merge_control':fingerprints[i]=key
                    elif fingerprints[i]!=key:raise ValueError('Unpaired oracle inputs')
                    amplitude=latent+delta*truth
                    threshold=(amplitude>delta/2).astype(int)
                    balanced=np.zeros(len(truth),int);balanced[np.argsort(amplitude,kind='stable')[-300:]]=1
                    raw_accuracy=accuracy(truth,labels,eligible);coverage=float(np.isin(labels,eligible).mean())
                    row=dict(stage=stage,replicate=i,groups=int(len(eligible)),raw_accuracy=raw_accuracy,
                        raw_coverage=coverage,raw_ARI=float(adjusted_rand_score(truth,labels)),
                        oracle_threshold_accuracy=float(np.mean(threshold==truth)),
                        oracle_balanced_accuracy=float(np.mean(balanced==truth)),
                        oracle_balanced_ARI=float(adjusted_rand_score(truth,balanced)),decisions={})
                    for method in ('raw','NB_calibrated_gate','state_calibrated_gate'):
                        matches=[r for r in evaluations if r['stage']==stage and r['condition']=='two_clones' and r['method']==method and r['replicate']==i]
                        if len(matches)!=1:raise ValueError('Missing locked decision')
                        accepted=bool(matches[0]['event'])
                        row['decisions'][method]=dict(accepted=accepted,
                            accuracy_all_cells_if_gate_applied=raw_accuracy if accepted else 0.,
                            coverage_if_gate_applied=coverage if accepted else 0.)
                    rows.append(row)
                if i+1!=200:raise ValueError('Incomplete alternative panel')
        summary=[]
        for stage in ('merge_control','exact_matched'):
            selected=[r for r in rows if r['stage']==stage]
            entry=dict(stage=stage,n_simulated_specimens=len(selected),
                raw_exactly_two=sum(r['groups']==2 for r in selected),metrics={})
            for metric in ('raw_accuracy','raw_coverage','raw_ARI','oracle_threshold_accuracy','oracle_balanced_accuracy','oracle_balanced_ARI'):
                values=np.asarray([r[metric] for r in selected])
                entry['metrics'][metric]=dict(mean=float(values.mean()),sd_between_simulated_specimens=float(values.std(ddof=1)),
                    percentile_2_5=float(np.quantile(values,.025)),percentile_97_5=float(np.quantile(values,.975)))
            entry['gates']={method:dict(accepted=sum(r['decisions'][method]['accepted'] for r in selected),
                mean_accuracy_all_cells_if_applied=float(np.mean([r['decisions'][method]['accuracy_all_cells_if_gate_applied'] for r in selected])),
                mean_coverage_if_applied=float(np.mean([r['decisions'][method]['coverage_if_gate_applied'] for r in selected])))
                for method in ('raw','NB_calibrated_gate','state_calibrated_gate')}
            summary.append(entry)
        (a.out/'records.json').write_text(json.dumps(rows,indent=2))
        (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,analytic_per_cell_Bayes_accuracy=analytic,
            independent_quadrature_accuracy=a0+a1,quadrature_error_estimate=e0+e1,
            interpretation='Analytic limit is expected accuracy with known equal priors and independent cell labels. Fixed-count top-half oracle is separate; empirical oracle results are not per-realization bounds.',
            replication='200 simulated specimens; cells are not independent validation specimens; displayed percentiles describe between-specimen spread, not confidence intervals',
            causal_limit='Bimodal transcriptional state can produce the same amplitude distribution; even oracle amplitude does not establish DNA causation',
            decision='No method/cutoff/criterion changes; no oracle input to inference'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
