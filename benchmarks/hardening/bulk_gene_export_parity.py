"""Compare public bulk output with its preserved wheel and audit gene effects."""
import argparse
import fcntl
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


CALL="""import sys
from haplodose.bulk import run_bulk
from pathlib import Path
data,out,mode=Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3]
run_bulk(data/'counts.tsv',data/'samples.tsv',data/'positions.tsv',out,
         n_programmes=2,bin_genes=10,bin_stat=mode,counts_kind='raw',genome_build='hg38',
         coordinate_system='zero-based-half-open',gene_namespace='custom',device='cpu')
"""


def contents(path):
    return gzip.open(path,'rb').read() if path.suffix=='.gz' else path.read_bytes()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    old_artifact=root/'dev/allele_cell_model_20261108/validation_23002524/installed_artifact.json'
    old=json.loads(old_artifact.read_text());old_python=Path(old['python'])
    # No concurrent installation or mutation of the preserved wheel environment.
    lock=(old_python.parent.parent/'.hardening-validate.lock').open('r')
    fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
    for name in ('wheel','native'):
        assert file_record(old[name]['path'])['sha256']==old[name]['sha256']
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
                  old_artifact=old_artifact,baseline=args.task/'baseline/source_manifest.json'),
                  dict(seed=20261108,modes=['sum','median']),
                  dict(scope='Public bulk regression and gene effect arithmetic; no biological validation'),
                  mode='bulk_gene_export_parity')
    try:
        data=args.out/'data';data.mkdir()
        rng=np.random.default_rng([20261108,1201]);genes=np.array([f'G{i:04d}' for i in range(120)])
        order=rng.permutation(120);names=['T2','001','N2','T0','N0','N1','T1']
        normal=np.array([False,True,True,False,True,True,False])
        base=rng.lognormal(3,.8,120);rates=np.tile(base[:,None],(1,7))
        rates[30:60,~normal]*=1.35;rates[90:,~normal]*=.65
        counts=pd.DataFrame(rng.poisson(rates*20)[order],index=genes[order],columns=names)
        counts.to_csv(data/'counts.tsv',sep='\t')
        pd.DataFrame(dict(sample=names,role=np.where(normal,'normal','tumour'))).to_csv(data/'samples.tsv',sep='\t',index=False)
        within=np.tile(np.arange(60),2)
        pd.DataFrame(dict(gene=genes,chromosome=np.repeat(['chr1','chr2'],60),
                          start=130000000+within*10000,end=130000500+within*10000)).to_csv(data/'positions.tsv',sep='\t',index=False)
        (args.out/'executed_call.py').write_text(CALL)
        records=[];new_effects=[]
        old_source=json.loads((args.task/'baseline/source_manifest.json').read_text())
        for mode in ('sum','median'):
            for label,python in (('old',old_python),('new',Path(sys.executable))):
                destination=args.out/(mode+'_'+label)
                result=subprocess.run([str(python),'-I','-c',CALL,str(data),str(destination),mode],
                                      capture_output=True,text=True,check=False,cwd=args.out)
                (args.out/(mode+'_'+label+'.log')).write_text(result.stdout+result.stderr)
                assert result.returncode==0,(mode,label,result.returncode)
                manifest=json.loads((destination/'run_manifest.json').read_text())
                assert manifest['status']=='complete' and manifest['backend']=='cpu'
                if label=='old':
                    for name,record in manifest['execution']['source'].items():
                        if name.endswith('.py'):
                            assert record['sha256']==old_source['python/haplodose/'+name]
                for record in manifest['output']['files'].values():
                    assert file_record(record['path'])['sha256']==record['sha256']
            previous=args.out/(mode+'_old');current=args.out/(mode+'_new')
            compared=[]
            for name in ('bins.tsv','genes.tsv.gz','segments.tsv','arms.tsv','samples.tsv','log_ratio.npy'):
                assert contents(previous/name)==contents(current/name),'Existing output changed: '+mode+'/'+name
                compared.append(name)
            before=json.loads((previous/'run.json').read_text());after=json.loads((current/'run.json').read_text())
            metadata=after.pop('gene_rna_log_ratio');assert before==after
            assert metadata['calibration_status']=='unvalidated' and metadata['programme_correction'] is False
            mapping=pd.read_csv(current/'genes.tsv.gz',sep='\t')
            effects=np.load(current/'gene_rna_log_ratio.npy');new_effects.append(effects)
            exported=pd.read_csv(current/'samples.tsv',sep='\t',dtype={'sample':str})
            provenance=pd.read_csv(current/'gene_ratio_metadata.tsv',sep='\t',dtype={'sample':str},float_precision='round_trip')
            assert exported['sample'].tolist()==names==provenance['sample'].tolist()
            assert effects.dtype==np.float64 and effects.shape==(7,120)
            groups=[np.flatnonzero(mapping.bin.to_numpy()==b) for b in sorted(mapping.bin.unique())]
            # Independent scalar sums and ratios, retaining the declared raw
            # measurement scale. The operation sequence matches the formula.
            for i,sample in enumerate(names):
                refs=[name for name,is_normal in zip(names,normal) if is_normal and name!=sample]
                total=sum(int(counts.loc[g,name]) for name in refs for g in mapping.gene)
                library=sum(int(counts.loc[g,sample]) for g in mapping.gene)
                effect=np.array([np.log((int(counts.loc[g,sample])+.5)/(library*(sum(int(counts.loc[g,name]) for name in refs)/total)+.5)) for g in mapping.gene])
                centre=float(np.median([np.median(effect[indices]) for indices in groups]))
                np.testing.assert_array_equal(effects[i],effect-centre)
                assert provenance.loc[i,'centre_rna_log_ratio']==centre
                assert provenance.loc[i,'reference_input_count']==total
                assert provenance.loc[i,'mapped_input_count']==library
                assert provenance.loc[i,'n_reference_samples']==len(refs)
            records.append(dict(mode=mode,old_outputs_identical=compared,run_metadata_unchanged_except_new_export=True,
                                scalar_gene_effects_checked=int(effects.size),axis_identity_exact=True))
        np.testing.assert_array_equal(*new_effects)
        write_json(dict(modes=records,raw_gene_effects_independent_of_bin_stat=True,
                        biological_validation=False,original_recovery_failures_resolved=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
