"""Independent canonical parity, unchanged criteria and conditional redraw summaries."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import beta,binom
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def data_parity(left,right):
    a,b=(json.loads((p/'summary.json').read_text()) for p in (left,right))
    assert len(a['cases'])==len(b['cases'])==5
    for x,y in zip(a['cases'],b['cases']):
        assert x['seed']==y['seed'] and x['criteria']==y['criteria']
        for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
            pd.testing.assert_frame_equal(pd.read_csv(left/f'seed_{x["seed"]}/fit'/name,sep='\t'),
                                          pd.read_csv(right/f'seed_{y["seed"]}/fit'/name,sep='\t'),check_exact=True)
        assert np.array_equal(np.load(left/f'seed_{x["seed"]}/fit/posteriors.npy'),
                              np.load(right/f'seed_{y["seed"]}/fit/posteriors.npy'),equal_nan=True)


def bounds(k,n):
    # Ten conditional method/case probabilities, two tails each.
    alpha=.05/20
    lower=0. if k==0 else float(beta.ppf(alpha,k,n-k+1))
    upper=1. if k==n else float(beta.ppf(1-alpha,k+1,n-k))
    lo=0. if k==0 else brentq(lambda p:binom.sf(k-1,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    hi=1. if k==n else brentq(lambda p:binom.cdf(k,n,p)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    assert abs(lo-lower)<=128*np.finfo(float).eps and abs(hi-upper)<=128*np.finfo(float).eps
    return [lower,upper]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--gate',action='store_true')
    p.add_argument('--gate-root',default='gate_v1')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    gate=a.task/a.gate_root
    gate_runs=json.loads((gate/'runs.json').read_text())
    files=[gate/'runs.json'] if a.gate else sorted(a.task.glob('chunk_*_v1/runs.json'))
    runs=[r for f in files for r in json.loads(f.read_text())]
    assert len(runs)==(4 if a.gate else 20)
    bank=json.loads((a.task/'seed_bank.json').read_text())
    for r in bank['draws']:
        value=int.from_bytes(hashlib.blake2b(f'HaploDose|split-draw-probe|20261109|{r["rep"]}'.encode(),digest_size=4).digest(),'little')
        assert value==r['seed'] and value!=bank['canonical_seed']
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',seeds=a.task/'seed_bank.json',
        **{str(i):f for i,f in enumerate(files)}),dict(gate=a.gate),
        dict(scope='Conditional count-draw diagnostic; biological n=0'),mode='split_draw_audit')
    argv=sys.argv
    try:
        for arm in ('bin','gene'):
            data_parity(gate/arm/'default/panel',gate/arm/'canonical/panel')
        data_parity(gate/'bin/default/panel',a.task.parent/'exact_profile_merge_20261109/production_v1/panel')
        # Recount the original six criteria independently, including every failure.
        helper=load('redraw_criteria_audit',Path(__file__).with_name('offset_recovery_audit.py'))
        arms=[str((Path(r['folder'])/'panel').relative_to(a.task)) for r in runs]
        sys.argv=[str(helper.__file__),'--out',str(a.out/'criteria'),'--task',str(a.task),'--arms',*arms]
        helper.main()
        canonical={};fold_count=0;gene_count=0
        for run in gate_runs:
            if not run['override_enabled']:continue
            for record in run['folds']:
                with np.load(record['path']) as z:
                    key=(int(z['case']),tuple(z['rows']))
                    if key in canonical:
                        assert np.array_equal(canonical[key]['full'],z['full']) and np.array_equal(canonical[key]['size'],z['size'])
                    canonical[key]=dict(full=z['full'],size=z['size'])
        assert len(canonical)==15
        for run in runs:
            expected=bank['canonical_seed'] if a.gate else bank['draws'][run['rep']]['seed']
            assert run['effective_seed']==expected
            summary=json.loads((Path(run['folder'])/'panel/summary.json').read_text())
            assert summary['params']['seed']==bank['canonical_seed'] and summary['stage']=='hard'
            lookup={}
            for record in run['folds']:
                with np.load(record['path']) as z:
                    key=(int(z['case']),tuple(z['rows']))
                    assert int(z['requested_seed'])==bank['canonical_seed'] and int(z['effective_seed'])==expected
                    assert float(z['share'])==.5 and np.isfinite(z['size']).all() and (z['size']>0).all()
                    assert np.array_equal(z['full'],canonical[key]['full']) and np.array_equal(z['size'],canonical[key]['size'])
                    assert np.array_equal(z['prop']+z['test'],z['full'])
                    assert (z['prop']>=0).all() and (z['prop']<=z['full']).all()
                    assert np.equal(z['prop'],np.floor(z['prop'])).all()
                    lookup[key]=dict(prop=z['prop'],test=z['test'])
                fold_count+=1
            if run['override_enabled']:assert len(lookup)==15
            for record in run['gene_folds']:
                with np.load(record['path']) as z:
                    key=(record['case_seed'],tuple(z['global_rows']))
                    assert int(z['seed'])==expected and float(z['gene_size'])==10.
                    if run['override_enabled']:
                        assert np.array_equal(z['prop_bins'],lookup[key]['prop'])
                        assert np.array_equal(z['full_bins']-z['prop_bins'],lookup[key]['test'])
                gene_count+=1
            # Genes, reference-derived initial dispersion and bin universe must
            # remain fixed when only the count-splitting seed changes.
            for case in summary['cases']:
                folder=Path(run['folder'])/f'panel/seed_{case["seed"]}/fit'
                previous=gate/f'{run["arm"]}/default/panel/seed_{case["seed"]}/fit'
                for name in ('genes.tsv.gz','bins.tsv'):
                    pd.testing.assert_frame_equal(pd.read_csv(folder/name,sep='\t'),pd.read_csv(previous/name,sep='\t'),check_exact=True)
                current=json.loads((folder/'run_manifest.json').read_text())
                old=json.loads((previous/'run_manifest.json').read_text())
                assert current['input']==old['input']
        conditional=[];paired=[]
        if not a.gate:
            assert {(r['arm'],r['rep']) for r in runs}=={(m,i) for m in ('bin','gene') for i in range(10)}
            scores={}
            for run in runs:
                summary=json.loads((Path(run['folder'])/'panel/summary.json').read_text())
                for row in summary['cases']:scores[(run['arm'],run['rep'],row['seed'])]=row
            cases=sorted({key[2] for key in scores})
            for case in cases:
                for arm in ('bin','gene'):
                    observations=[scores[(arm,i,case)] for i in range(10)]
                    k=sum(r['original_assertion_passed'] for r in observations)
                    conditional.append(dict(case=case,arm=arm,passed=k,draws=10,
                        conditional_draw_probability_interval=bounds(k,10),
                        pass_by_draw=[r['original_assertion_passed'] for r in observations],
                        failed_criteria_by_draw=[[name for name,value in r['criteria'].items() if not value['passed']] for r in observations]))
                paired.append(dict(case=case,bin_only=sum(scores[('bin',i,case)]['original_assertion_passed'] and not scores[('gene',i,case)]['original_assertion_passed'] for i in range(10)),
                    gene_only=sum(scores[('gene',i,case)]['original_assertion_passed'] and not scores[('bin',i,case)]['original_assertion_passed'] for i in range(10))))
        write_json(dict(all_checks_passed=True,runs=len(runs),original_criteria_checked=30*len(runs),
            fold_models_checked=fold_count,gene_fold_models_checked=gene_count,conditional=conditional,paired=paired,
            biological_specimens=0,production_changed=False,scientific_criteria_changed=False,
            interpretation='Technical redraw probabilities conditional on5 exposed matrices, not biological error control or a selected best-seed procedure'),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU independent original criteria, conserved folds, exact features/inputs and binomial inversion')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
