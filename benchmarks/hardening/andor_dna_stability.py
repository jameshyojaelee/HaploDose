"""Complete DNA refits on prespecified complementary read-count splits."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--start',type=int,required=True)
    parser.add_argument('--stop',type=int,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    assert 0<=args.start<args.stop<=20
    helper=Path(__file__).with_name('andor_dna_stability_common.py')
    spec=importlib.util.spec_from_file_location('dna_common',helper)
    common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
    counts=args.task/'count_audit_v1';reconstruction=args.task/'reconstruction_v1'
    common.verify_run(counts);common.verify_run(reconstruction)
    module=common.dna_module()
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,DNA_source=Path(module.__file__),
                  design=args.task/'DESIGN.md',counts=counts/'counts.npz',
                  annotation=counts/'annotation.tsv',bin_order=counts/'bins.tsv',header=counts/'read_groups.txt',
                  count_gate=counts/'run_manifest.json',reconstruction=reconstruction/'run_manifest.json'),
                  dict(start=args.start,stop=args.stop,allocation=.5,master_seed=20261108),
                  dict(scope='Conditional sequencing sensitivity in one exposed DNA specimen'),mode='andor_dna_stability')
    try:
        with np.load(counts/'counts.npz') as z:
            full=z['new'];cells=z['cells']
        bins=pd.read_csv(counts/'bins.tsv',sep='\t');results=[]
        for rep in range(args.start,args.stop):
            pair=args.out/f'pair_{rep:02d}';pair.mkdir()
            rng=np.random.default_rng(np.random.SeedSequence([20261108,1401,rep]))
            first=rng.binomial(full,.5);second=full-first
            assert np.all(first>=0) and np.all(second>=0)
            np.testing.assert_array_equal(first+second,full)
            fits=[];inputs=[]
            for label,matrix in (('A',first),('B',second)):
                parent=pair/label;parent.mkdir();folder=parent/'bins';folder.mkdir()
                shutil.copy2(counts/'annotation.tsv',parent/'bin_annotation_500kb.tsv')
                shutil.copy2(counts/'read_groups.txt',folder/'SNU-601.rg.txt')
                common.write_counts(folder/'SNU-601.tsv.gz',matrix,cells,bins)
                inputs.append(file_record(folder/'SNU-601.tsv.gz'))
                fits.append(common.fit(module,folder,parent/'DNA'))
            comparisons=dict(A_to_full=common.compare(reconstruction/'corrected/SNU-601',fits[0]),
                             B_to_full=common.compare(reconstruction/'corrected/SNU-601',fits[1]),
                             A_to_B=common.compare(fits[0],fits[1]))
            result=dict(pair=rep,seed=[20261108,1401,rep],half_counts=inputs,
                        full_fragments=int(full.sum()),A_fragments=int(first.sum()),B_fragments=int(second.sum()),
                        comparisons=comparisons)
            write_json(result,pair/'comparison.json');results.append(result)
            print(json.dumps(dict(pair=rep,status='complete')),flush=True)
        write_json(dict(pairs=results,biological_specimens=1,
                        uncertainty='Conditional sequencing perturbations, complementary halves dependent',
                        no_RNA_input=True,validation_claim=False),args.out/'report.json')
        common.write_inventory(args.out)
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
