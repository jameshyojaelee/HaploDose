"""Fixed-membership perturbation test for SNP confirmation entering discovery."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import jax
import numpy as np
import pandas as pd

from haplodose.alleles import AlleleParams, SpecimenAlleles, balance_loglik_ratio, state_b
from haplodose.clones import CloneParams, SpecimenModel
from haplodose.provenance import start_run, finish_run


def digest(array):
    array=np.ascontiguousarray(array)
    return hashlib.sha256(str(array.dtype).encode()+json.dumps(array.shape).encode()+array.tobytes()).hexdigest()


def run_case(a,d,snps,proposal_only):
    index=np.arange(0,40,2) if proposal_only else np.arange(40)
    allele=SpecimenAlleles(a[:,index],d[:,index],snps.iloc[index],50.,AlleleParams())
    pools=allele.pool(np.ones((120,1)))
    if not proposal_only:
        model=SimpleNamespace(alleles=allele,params=CloneParams())
        SpecimenModel.set_allele_pools(model,*pools)
        proposals=SpecimenModel.allele_splits(model,0,np.arange(120))
        post=model.allele_post[0]
    else:
        # Physical subset avoids counting zero-depth held-out SNPs toward
        # the unchanged minimum-SNP rule. Retain original genomic distances.
        post=allele.fit(*pools)[0][0]
        proposals=[]
        for r in allele.runs(post,pools[1][0]):
            if not r['resolved']:
                continue
            selected=np.arange(r['first_snp'],r['last_snp']+1)
            b=.5+allele.params.cell_shrink*(post[selected]@state_b()-.5)
            llr=balance_loglik_ratio(allele.A[:,selected],allele.D[:,selected],b,allele.params.cell_kappa)
            labels=(llr>0).astype(int)
            if min(labels.sum(),len(labels)-labels.sum())>=CloneParams().min_cells:
                proposals.append((labels,f"alleles_{r['chromosome']}:{r['start']}-{r['end']}"))
    even_post=post if proposal_only else post[::2]
    return dict(nu=float(allele.nu),kappa=float(allele.kappa),
                proposal_posterior_sha256=digest(even_post),
                runs=allele.runs(post,pools[1][0]),
                proposals=[dict(labels=g.tolist(),source=src) for g,src in proposals]),even_post


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    assert jax.default_backend()=='cpu'
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),
                  dict(cells=120,snps=40,reference_kappa=50.,odd_counts=[10,19,1]),
                  dict(scope='Fixed memberships and reference; local leakage test, not independent workflow validation'),
                  mode='allele_confirmation_leakage')
    try:
        snps=pd.DataFrame(dict(chromosome='chr1',pos=np.arange(40)*1000+1,cM=np.arange(40)*.001))
        d=np.full((120,40),20.)
        a=np.repeat(np.repeat([10.,12.,16.],40)[:,None],40,axis=1)
        results={};posteriors={};count_hashes=[]
        for name,count in (('odd_balanced',10.),('odd_same_orientation',19.),('odd_reverse_orientation',1.)):
            changed=a.copy();changed[:,1::2]=count
            count_hashes.append(dict(case=name,proposal_A=digest(changed[:,::2]),
                                     proposal_D=digest(d[:,::2]),confirmation_A=digest(changed[:,1::2])))
            np.savez_compressed(args.out/(name+'_inputs.npz'),A=changed,D=d)
            for arm in ('current','proposal_only_reference'):
                record,post=run_case(changed,d,snps,arm!='current')
                results[name+':'+arm]=record
                posteriors[name+':'+arm]=post
        assert len({r['proposal_A'] for r in count_hashes})==1
        assert len({r['proposal_D'] for r in count_hashes})==1
        assert len({r['confirmation_A'] for r in count_hashes})==3
        names=[r['case'] for r in count_hashes]
        comparisons={}
        for arm in ('current','proposal_only_reference'):
            baseline=results[names[0]+':'+arm]
            comparisons[arm]=[]
            for name in names[1:]:
                record=results[name+':'+arm]
                difference=float(np.max(abs(posteriors[name+':'+arm]-posteriors[names[0]+':'+arm])))
                comparisons[arm].append(dict(contrast=name+' minus '+names[0],
                    posterior_max_difference=difference,
                    fitted_parameters_identical=(record['nu'],record['kappa'])==(baseline['nu'],baseline['kappa']),
                    proposal_posteriors_identical=record['proposal_posterior_sha256']==baseline['proposal_posterior_sha256'],
                    runs_identical=record['runs']==baseline['runs'],
                    proposals_identical=record['proposals']==baseline['proposals']))
        assert all(r['proposal_posteriors_identical'] and r['fitted_parameters_identical']
                   and r['runs_identical'] and r['proposals_identical']
                   for r in comparisons['proposal_only_reference']), 'Held-out data entered physically subset reference'
        np.savez(args.out/'posteriors.npz',**posteriors)
        (args.out/'details.json').write_text(json.dumps(results,indent=2))
        (args.out/'report.json').write_text(json.dumps(dict(inputs=count_hashes,comparisons=comparisons,
             biological_specimens=0,production_changed=False,
             scope='Local allele proposal at fixed memberships; upstream memberships and adaptive selection not audited'),indent=2))
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
