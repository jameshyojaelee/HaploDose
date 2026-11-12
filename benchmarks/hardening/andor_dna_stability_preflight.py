"""Exercise DNA refits, count identity and label invariance before real outcomes."""
import argparse
from collections import Counter
import gzip
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def sibling(name):
    path=Path(__file__).with_name(name+'.py')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    common=sibling('andor_dna_stability_common');counter=sibling('andor_count_repair_audit')
    module=common.dna_module()
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
                  helper=Path(common.__file__),counter=Path(counter.__file__),DNA_source=Path(module.__file__)),
                  dict(seed=20261108,cells=48,bins=240),dict(scope='Synthetic DNA computation, biological n0'),mode='andor_dna_preflight')
    try:
        np.testing.assert_array_equal(counter.boundary_flow(np.array([[2,3,2]]),np.array([[3,3,1]])),[[0,1,1]])
        rejected=0
        for old,new in (([[1,0]],[[0,1]]),([[1,0]],[[1,1]]),([[2,0,0]],[[1,1,0]])):
            try:
                counter.boundary_flow(np.array(old),np.array(new))
            except AssertionError:
                rejected+=1
            else:
                raise AssertionError('Invalid correction passed boundary-flow gate')
        reader_bins=pd.DataFrame(dict(chromosome=['chr1','chr1','chr1','chr2'],bin=[0,1,2,0]))
        reader_cells=np.array(['CELL-B','CELL-A']);reader_matrix=np.array([[0,3,0,7],[8,0,2,0]])
        reader_path=args.out/'reader_fixture.tsv.gz'
        common.write_counts(reader_path,reader_matrix,reader_cells,reader_bins)
        with gzip.open(reader_path,'at') as stream:
            stream.write('CELL-A\tchrX\t3\t9\n')
        reader_observed,reader_totals,_=counter.read_counts(reader_path,reader_cells,reader_bins)
        np.testing.assert_array_equal(reader_observed,reader_matrix)
        assert reader_totals==Counter({('CELL-B','chr1'):3,('CELL-B','chr2'):7,
                                     ('CELL-A','chr1'):10,('CELL-A','chrX'):9})
        rng=np.random.default_rng([20261108,1400]);cells=np.array([f'CELL-{i:03d}' for i in range(48)])
        within=np.tile(np.arange(120),2)
        ann=pd.DataFrame(dict(chromosome=np.repeat(['chr1','chr2'],120),bin=260+within,
                             start=130000000+500000*within,end=130500000+500000*within,
                             arm=np.repeat(['chr1q','chr2q'],120),gc=rng.uniform(.3,.7,240),
                             n_frac=0.,blacklist=0.,mappability=1.))
        rate=np.broadcast_to(np.exp(2*(ann.gc.to_numpy()-.5)),(48,240)).copy()
        rate[:24,:30]*=1.5;rate[24:,120:150]*=.5
        full=rng.poisson(rate/rate.sum(axis=1,keepdims=True)*300000)
        a=rng.binomial(full,.5);b=full-a;np.testing.assert_array_equal(a+b,full)
        fits=[]
        for label,matrix in (('full',full),('A',a),('B',b)):
            parent=args.out/label;parent.mkdir();folder=parent/'bins';folder.mkdir()
            ann.to_csv(parent/'bin_annotation_500kb.tsv',sep='\t',index=False)
            (folder/'SNU-601.rg.txt').write_text(''.join(f'@RG\tID:{c}\tSM:{c}\n' for c in cells))
            common.write_counts(folder/'SNU-601.tsv.gz',matrix,cells,ann[['chromosome','bin']])
            loaded_ids,loaded=module.load_counts(folder/'SNU-601.tsv.gz',ann)
            np.testing.assert_array_equal(loaded_ids,cells);np.testing.assert_array_equal(loaded,matrix)
            fits.append(common.fit(module,folder,parent/'DNA'))
        self_score=common.compare(fits[0],fits[0]);assert self_score['mean_absolute_relative_profile_change']==0.
        assert self_score['fraction_integer_CN_entries_changed']==0.
        assert self_score['one_to_one_matched_assigned']==self_score['anchor_assigned']
        renamed=args.out/'renamed';shutil.copytree(fits[0],renamed)
        frame=pd.read_csv(renamed/'cells.tsv',sep='\t');qc=json.loads((renamed/'qc.json').read_text())
        names=list(qc['clones']);rename={name:'PERMUTED_'+str(len(names)-i) for i,name in enumerate(names)}
        frame['clone']=frame['clone'].replace(rename)
        frame.iloc[rng.permutation(len(frame))].to_csv(renamed/'cells.tsv',sep='\t',index=False)
        qc['clones']={rename[k]:v for k,v in qc['clones'].items()};write_json(qc,renamed/'qc.json')
        with np.load(renamed/'cn.npz',allow_pickle=True) as z:
            values={k:z[k] for k in z.files}
        permutation=rng.permutation(len(values['cells']))
        for name in ('cells','cn','seg_ratio'):
            values[name]=values[name][permutation]
        np.savez_compressed(renamed/'cn.npz',**values)
        for filename in ('clone_bins.tsv','clone_arms.tsv','clone_differences.tsv'):
            table=pd.read_csv(renamed/filename,sep='\t').rename(columns=rename)
            for name in ('clone_a','clone_b'):
                if name in table:
                    table[name]=table[name].replace(rename)
            table.to_csv(renamed/filename,sep='\t',index=False)
        assert common.compare(fits[0],renamed)==self_score,'DNA agreement depends on label or row order'
        write_json(dict(negative_flow_controls_rejected=rejected,complete_DNA_fits=3,
                        independent_integer_reader_identity_and_zero_fill=True,
                        exact_count_reader=True,label_and_row_order_invariant=True,
                        self_comparison=self_score,half_comparison=common.compare(fits[1],fits[2]),
                        biological_specimens=0),args.out/'report.json')
        common.write_inventory(args.out)
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
