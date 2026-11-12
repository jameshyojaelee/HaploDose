"""Descriptive DNA coassignment stability without choosing a label matching or cutoff."""
import argparse
from collections import Counter
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record

UNRESOLVED={'excluded','unassigned'}


def pair_counts(anchor,candidate):
    """Fixed-reference cell-pair counts; unresolved labels never count as a group."""
    assert len(anchor)==len(candidate)
    names=sorted(set(anchor)-UNRESOLVED)
    groups={name:Counter(right for left,right in zip(anchor,candidate)
                         if left==name and right not in UNRESOLVED) for name in names}
    sizes=Counter(anchor);result=[]
    for left,right in itertools.combinations_with_replacement(names,2):
        a,b=groups[left],groups[right]
        if left==right:
            total=sizes[left]*(sizes[left]-1)//2
            assigned=sum(a.values());both=assigned*(assigned-1)//2
            coassigned=sum(n*(n-1)//2 for n in a.values())
        else:
            total=sizes[left]*sizes[right];both=sum(a.values())*sum(b.values())
            coassigned=sum(n*b.get(key,0) for key,n in a.items())
        assert 0<=coassigned<=both<=total
        result.append(dict(anchor_left=left,anchor_right=right,within_reference_group=left==right,
                           fixed_reference_pairs=total,both_assigned_pairs=both,
                           coassigned_pairs=coassigned,separated_pairs=both-coassigned,
                           unresolved_pairs=total-both))
    return result


def read_groups(path):
    frame=pd.read_csv(path,sep='\t',dtype={'rg':str,'clone':str},keep_default_na=False)
    assert (frame.rg!='').all() and frame.rg.is_unique and (frame.clone!='').all()
    return frame.set_index('rg').clone


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('andor_dna_stability_audit.py')
    spec=importlib.util.spec_from_file_location('verified_dna',helper)
    verified=importlib.util.module_from_spec(spec);spec.loader.exec_module(verified)
    inputs=dict(script=Path(__file__),helper=helper,design=args.task/'CONTRAST_DESIGN.md',
                preflight=args.task/'contrast_preflight_v1/run_manifest.json',
                audit=args.task/'audit_v1/run_manifest.json',
                anchor=args.task/'reconstruction_v1/corrected/SNU-601/cells.tsv')
    run=start_run(args.out,inputs,{},dict(scope='Exploratory conditional stability, one exposed DNA specimen'),
                  mode='andor_dna_contrast_stability')
    try:
        verified.verified_run(args.task/'contrast_preflight_v1')
        verified.verified_run(args.task/'audit_v1')
        verified.verified_run(args.task/'reconstruction_v1')
        anchor=read_groups(inputs['anchor']);rows=[];coverage=[]
        for worker in range(4):
            directory=args.task/f'worker_{worker}_v1';verified.verified_run(directory)
            for rep in range(5*worker,5*(worker+1)):
                for half in ('A','B'):
                    path=directory/f'pair_{rep:02d}/{half}/DNA/SNU-601/cells.tsv'
                    candidate=read_groups(path)
                    assert set(anchor.index)==set(candidate.index)
                    candidate=candidate.loc[anchor.index]
                    records=pair_counts(anchor.tolist(),candidate.tolist())
                    for record in records:
                        record.update(pair=rep,half=half)
                    rows.extend(records)
                    for name in sorted(set(anchor)-UNRESOLVED):
                        group=candidate.loc[anchor==name];counts=Counter(group)
                        coverage.append(dict(pair=rep,half=half,anchor_group=name,reference_cells=len(group),
                                             assigned_cells=sum(value for key,value in counts.items() if key not in UNRESOLVED),
                                             excluded_cells=counts['excluded'],unassigned_cells=counts['unassigned']))
        table=pd.DataFrame(rows)
        table.to_csv(args.out/'cell_pair_counts.tsv',sep='\t',index=False)
        pd.DataFrame(coverage).to_csv(args.out/'group_coverage.tsv',sep='\t',index=False)
        summaries=[]
        for (left,right),group in table.groupby(['anchor_left','anchor_right'],sort=True):
            assert len(group)==40 and group.fixed_reference_pairs.nunique()==1
            total=int(group.fixed_reference_pairs.iloc[0]);assert total>0
            # Each half is a conditional perturbation; no binomial confidence claim.
            coassigned=group.coassigned_pairs.to_numpy()/total
            separated=group.separated_pairs.to_numpy()/total
            unresolved=group.unresolved_pairs.to_numpy()/total
            assert (group.coassigned_pairs+group.separated_pairs+group.unresolved_pairs==group.fixed_reference_pairs).all()
            summaries.append(dict(anchor_left=left,anchor_right=right,within_reference_group=left==right,
                                  fixed_reference_pairs=total,conditional_half_fits=40,independent_pair_seeds=20,
                                  mean_coassigned_fraction=float(np.mean(coassigned)),
                                  min_coassigned_fraction=float(np.min(coassigned)),max_coassigned_fraction=float(np.max(coassigned)),
                                  mean_separated_fraction=float(np.mean(separated)),
                                  mean_unresolved_fraction=float(np.mean(unresolved))))
        pd.DataFrame(summaries).to_csv(args.out/'contrast_summary.tsv',sep='\t',index=False)
        write_json(dict(reference_groups=len(set(anchor)-UNRESOLVED),reference_assigned_cells=int((~anchor.isin(UNRESOLVED)).sum()),
                        assay_cells=len(anchor),half_fits=40,independent_pair_seeds=20,
                        biological_specimens=1,contrasts=len(summaries),
                        matching='None; counts invariant to candidate group relabeling',
                        estimand='Conditional coassignment/separation of fixed full-depth reference cell pairs',
                        limitation='Full-depth labels use the same molecules; not independent confirmation or certified true clones',
                        old_arm_summary='Any inferred clone pair anywhere on an arm; not stability of a particular distinction',
                        new_support_cutoff=None,new_confidence_claim=False),args.out/'report.json')
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path.name not in ('run_manifest.json','output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
