"""Independently audit coordinate-repair impact before read-split results exist."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def load(name):
    path=Path(__file__).with_name(name+'.py')
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module,path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    audit,helper=load('andor_dna_stability_audit')
    metric,reference=load('andor_dna_independent')
    counts=args.task/'count_audit_v1';fit=args.task/'reconstruction_v1'
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,reference=reference,
                  count_manifest=counts/'run_manifest.json',fit_manifest=fit/'run_manifest.json',
                  design=args.task/'AUDIT_DESIGN.md'),{},dict(scope='One development specimen, coordinate-repair mechanism'),
                  mode='andor_dna_reconstruction_audit')
    try:
        manifest=audit.verified_run(counts);audit.verified_run(fit)
        with np.load(counts/'counts.npz') as archive:
            saved={key:archive[key] for key in archive.files}
        cells=saved['cells'];bins=pd.read_csv(counts/'bins.tsv',sep='\t')
        old,ot=audit.independent_counts(Path(manifest['input']['files']['old']['path']),cells,bins)
        new,nt=audit.independent_counts(Path(manifest['input']['files']['new']['path']),cells,bins)
        assert ot==nt
        np.testing.assert_array_equal(old,saved['old']);np.testing.assert_array_equal(new,saved['new'])
        ann_path=Path(manifest['input']['files']['annotation']['path'])
        assert ann_path.read_bytes()==(counts/'annotation.tsv').read_bytes()
        ann=pd.read_csv(ann_path,sep='\t')
        flow=audit.callable_flow(old,new,saved['flow'],bins,ann)
        comparison=json.loads((fit/'comparison.json').read_text())
        check=metric.audit_comparison(fit/'old/SNU-601',fit/'corrected/SNU-601',comparison)
        before=pd.read_csv(fit/'old/SNU-601/cells.tsv',sep='\t',dtype={'rg':str}).set_index('rg')
        after=pd.read_csv(fit/'corrected/SNU-601/cells.tsv',sep='\t',dtype={'rg':str}).set_index('rg')
        assert set(before.index)==set(after.index)==set(cells)
        pd.crosstab(before.clone,after.clone).to_csv(args.out/'group_transitions.tsv',sep='\t')
        callable_keys=set(zip(ann.loc[(ann.n_frac<=.1)&(ann.blacklist<=.1)&(ann.mappability>=.5),'chromosome'],
                              ann.loc[(ann.n_frac<=.1)&(ann.blacklist<=.1)&(ann.mappability>=.5),'bin']))
        keep=np.array([(row.chromosome,row.bin) in callable_keys for row in bins.itertuples()])
        changed=np.any(old[:,keep]!=new[:,keep],axis=1)
        trace=pd.DataFrame(dict(rg=cells,callable_counts_changed=changed,
                               boundary_fragments=saved['flow'].sum(axis=1),**flow)).set_index('rg')
        trace['old_group']=before.clone;trace['corrected_group']=after.clone
        # Only trusted locally produced archives; independently audited IDs are checked above.
        with np.load(fit/'old/SNU-601/cn.npz',allow_pickle=True) as z:
            a={key:z[key] for key in ('cells','bins','seg_ratio','cn')}
        with np.load(fit/'corrected/SNU-601/cn.npz',allow_pickle=True) as z:
            b={key:z[key] for key in ('cells','bins','seg_ratio','cn')}
        np.testing.assert_array_equal(a['bins'],b['bins'])
        ai={cell:i for i,cell in enumerate(a['cells'])};bi={cell:i for i,cell in enumerate(b['cells'])}
        unchanged_checks=0
        for cell in sorted(set(ai)&set(bi)):
            x=a['seg_ratio'][ai[cell]];y=b['seg_ratio'][bi[cell]]
            trace.loc[cell,'profile_entries_changed']=int(np.count_nonzero(x!=y))
            trace.loc[cell,'maximum_absolute_relative_profile_change']=float(np.max(np.abs(x.astype(float)-y.astype(float))))
            trace.loc[cell,'integer_CN_entries_changed']=int(np.count_nonzero(a['cn'][ai[cell]]!=b['cn'][bi[cell]]))
            if not trace.loc[cell,'callable_counts_changed']:
                np.testing.assert_array_equal(x,y);unchanged_checks+=1
        trace.to_csv(args.out/'cell_trace.tsv',sep='\t')
        write_json(dict(independent_numerical_check=check,count_matrices_independently_reconstructed=True,
                        callable_flow_totals={key:int(value.sum()) for key,value in flow.items()},
                        cells_with_callable_count_changes=int(changed.sum()),
                        unchanged_callable_count_cells_with_identical_profiles=unchanged_checks,
                        old_corrected_comparison=comparison,biological_specimens=1,
                        interpretation='DNA measurement sensitivity; no true-clone or RNA membership validation'),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
