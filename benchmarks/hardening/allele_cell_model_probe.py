"""Independent PMF counterexamples for allele proposal/selection model consistency."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.special import expit

from haplodose.alleles import AlleleParams, SpecimenAlleles, balanced_cells
from haplodose.clones import CloneParams, SpecimenModel
from haplodose.provenance import start_run, finish_run


def check(kappa):
    a = np.array([[19,0,800,0],[1,0,1400,0],[10,0,1000,0]],dtype=float)
    d = np.array([[20,0,2000,0]]*3,dtype=float)
    snps = pd.DataFrame(dict(chromosome='chr1',pos=np.arange(1,5),cM=np.arange(4)/100))
    params = AlleleParams(cell_kappa=kappa)
    allele = SpecimenAlleles(a,d,snps,50.,params)
    post = np.zeros((4,9));post[:,5] = 1.
    allele.runs = lambda q,n: [dict(resolved=True,first_snp=0,last_snp=3,chromosome='chr1',start=0,end=4)]
    model = SimpleNamespace(alleles=allele,allele_post=post[None],Dk=d.sum(axis=0)[None],
                            params=CloneParams(min_cells=1))
    proposals = SpecimenModel.allele_splits(model,0,np.arange(3))
    selected = proposals[0][0].tolist() if proposals else None
    def contrast(aa,dd):
        if kappa is None:
            return (stats.binom.logpmf(aa,dd,.77)-stats.binom.logpmf(aa,dd,.5)).sum(axis=1)
        return (stats.betabinom.logpmf(aa,dd,kappa*.77,kappa*.23)
                -stats.betabinom.logpmf(aa,dd,kappa*.5,kappa*.5)).sum(axis=1)
    llr = contrast(a[:,::2],d[:,::2])
    expected = (llr>0).astype(int).tolist()
    # Independently held fixed pooled fractions/run isolate the observation
    # model, without relying on an HMM fit or claiming a calibrated interval.
    balance = SpecimenAlleles(sparse.csr_matrix([[500.]]),sparse.csr_matrix([[1000.]]),snps.iloc[:1],50.,params)
    balance.fit = lambda ak,dk: (post[None,:1],np.array([[.8]]))
    balance.runs = lambda q,n: [dict(resolved=True,theta=.3,first_snp=0,last_snp=0)]
    got = balanced_cells(balance)[0].tolist()
    p_balanced = expit(-contrast(np.array([[500.]]),np.array([[1000.]]))[0])
    expected_balanced = [0] if p_balanced>=.99 else []
    return dict(kappa=kappa,proposal=selected,expected_proposal=expected,
                independent_log_odds=llr.tolist(),balanced=got,expected_balanced=expected_balanced,
                independent_balance_probability=float(p_balanced),
                proposal_matches=selected==expected,balance_matches=got==expected_balanced)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--expect',choices=('old','corrected'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    run=start_run(args.out,dict(script=Path(__file__)),dict(expect=args.expect),
                  dict(scope='Fixed independent PMF counterexample, not scientific calibration'),mode='allele_cell_model_probe')
    try:
        result=[check(None),check(2.)]
        assert result[0]['proposal_matches'] and result[0]['balance_matches'], 'Legacy fixture failed'
        correct=all(r['proposal_matches'] and r['balance_matches'] for r in result)
        (args.out/'report.json').write_text(json.dumps(dict(checks=result,model_consistent=correct),indent=2))
        assert correct==(args.expect=='corrected'), 'Observation-model contract differs from expected version'
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
