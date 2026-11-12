"""Separate PDF panels of already-defined conditional DNA stability measurements."""
import argparse
from collections import Counter
import importlib.util
import itertools
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def render(metrics,pairs,destination):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'font.weight':'normal',
        'axes.titlesize':6,'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,
        'legend.fontsize':6,'legend.frameon':False,'axes.spines.top':False,
        'axes.spines.right':False,'axes.linewidth':.5,'pdf.fonttype':42,'ps.fonttype':42})
    metrics=metrics.copy();metrics['half']=metrics.comparison.str[0]
    assert set(metrics.half)=={'A','B'}
    assert metrics.fixed_full_DNA_cells.nunique()==1
    total=int(metrics.fixed_full_DNA_cells.iloc[0])
    assert (metrics.one_to_one_matched_assigned<=metrics.candidate_assigned).all()
    assert (metrics.candidate_assigned<=total).all()
    assert not metrics.duplicated(['pair','half']).any()
    markers={'A':'o','B':'^'}

    def save(fig,name):
        fig.tight_layout();fig.savefig(destination/name,format='pdf');plt.close(fig)

    # KEY MESSAGE: Full-reference agreement must retain the assay coverage denominator.
    fig,ax=plt.subplots(figsize=(3.5,3.1));ax.plot([0,1],[0,1],color='#9E9E9E',lw=.6)
    for _,group in metrics.groupby('pair'):
        assert set(group.half)=={'A','B'}
        ax.plot(group.candidate_assigned/total,group.one_to_one_matched_assigned/total,color='#9E9E9E',lw=.5)
    for half,group in metrics.groupby('half'):
        ax.scatter(group.candidate_assigned/total,group.one_to_one_matched_assigned/total,
                   marker=markers[half],s=12,c='black',label=f'Half {half}')
    ax.set(xlim=(0,1.02),ylim=(0,1.02),xlabel='Assigned DNA cells / all full-depth assay cells',
           ylabel='Matched assignments / all full-depth assay cells',title='One DNA specimen; full-depth groups are a reference')
    ax.legend(loc='upper left');save(fig,'assignment_agreement.pdf')

    # KEY MESSAGE: Refit integer-CN scales are measurements, not known absolute ploidy.
    fig,ax=plt.subplots(figsize=(3.5,2.6))
    assert metrics.anchor_ploidy_scale.nunique()==1
    ax.axhline(metrics.anchor_ploidy_scale.iloc[0],color='#9E9E9E',lw=.8,label='Full-depth reference')
    for _,group in metrics.groupby('pair'):
        ax.plot(group.pair,group.candidate_ploidy_scale,color='#9E9E9E',lw=.5)
    for half,group in metrics.groupby('half'):
        ax.scatter(group.pair,group.candidate_ploidy_scale,marker=markers[half],s=12,c='black',label=f'Half {half}')
    ax.set(xlabel='Read-split seed index',ylabel='Fitted integer-CN scale (not absolute ploidy)',title='Complementary halves share one specimen')
    ax.legend(loc='best');save(fig,'integer_CN_scale.pdf')

    # Reference group numbers encode descending full-depth size; retain this order in all panels.
    names=sorted(set(pairs.anchor_left)|set(pairs.anchor_right),key=lambda name:int(name.rsplit('_d',1)[1]))
    within=pairs.loc[pairs.within_reference_group].copy()
    assert (within.fixed_reference_pairs>0).all()
    # KEY MESSAGE: Within-reference coassignment retains unresolved and separated pairs in its denominator.
    fig,ax=plt.subplots(figsize=(3.5,4.1));seeds=sorted(set(within.pair))
    for index,name in enumerate(names):
        group=within.loc[within.anchor_left==name]
        for half,half_group in group.groupby('half'):
            offsets=[(seeds.index(rep)-(len(seeds)-1)/2)*.6/max(len(seeds),1) for rep in half_group.pair]
            ax.scatter(half_group.coassigned_pairs/half_group.fixed_reference_pairs,index+np.asarray(offsets),
                       marker=markers[half],s=7,c='black',alpha=.55)
    ax.set(xlim=(0,1.02),yticks=range(len(names)),yticklabels=[name.rsplit('_',1)[1] for name in names],
           xlabel='Coassigned pairs / fixed within-reference cell pairs',ylabel='Full-depth inferred DNA group',
           title='All split halves; unresolved pairs remain in denominator')
    ax.invert_yaxis();save(fig,'within_group_coassignment.pdf')

    # KEY MESSAGE: Between-reference coassignment and unresolved pairs must be inspected separately.
    for field,title,filename in (
        ('coassigned_pairs','Between-reference coassignment / fixed cell pairs','between_group_coassignment.pdf'),
        ('unresolved_pairs','Unresolved pairs / fixed reference cell pairs','between_group_unresolved.pdf')):
        matrix=np.full((len(names),len(names)),np.nan)
        for (left,right),group in pairs.loc[~pairs.within_reference_group].groupby(['anchor_left','anchor_right']):
            assert group.fixed_reference_pairs.nunique()==1 and (group.fixed_reference_pairs>0).all()
            value=float(np.mean(group[field]/group.fixed_reference_pairs));assert 0<=value<=1
            i,j=names.index(left),names.index(right);matrix[min(i,j),max(i,j)]=value
        # Fixed [0,1] probability scale preserves zeros/ones across the two panels; no percentile clipping.
        fig,ax=plt.subplots(figsize=(3.5,3.2));shown=ax.imshow(matrix,vmin=0,vmax=1,cmap='Greys',interpolation='nearest')
        labels=[name.rsplit('_',1)[1] for name in names]
        ax.set(xticks=range(len(names)),yticks=range(len(names)),xticklabels=labels,yticklabels=labels,
               xlabel='Full-depth inferred DNA group',ylabel='Full-depth inferred DNA group',title=title)
        ax.tick_params(axis='x',rotation=90);fig.colorbar(shown,ax=ax,fraction=.045,pad=.03,label='Mean across half-depth fits')
        save(fig,filename)


def fixture():
    metrics=[];pairs=[];anchor=['SNU-601_d1']*3+['SNU-601_d2']*3
    for rep,half in itertools.product(range(2),('A','B')):
        candidate=['a','a','b','b','b','unassigned'];counter=Counter()
        for i,j in itertools.combinations(range(6),2):
            left,right=sorted((anchor[i],anchor[j]));valid=candidate[i]!='unassigned' and candidate[j]!='unassigned'
            counter[left,right,'total']+=1
            counter[left,right,'coassigned']+=int(valid and candidate[i]==candidate[j])
            counter[left,right,'unresolved']+=int(not valid)
        for left,right in itertools.combinations_with_replacement(sorted(set(anchor)),2):
            pairs.append(dict(pair=rep,half=half,anchor_left=left,anchor_right=right,within_reference_group=left==right,
                fixed_reference_pairs=counter[left,right,'total'],coassigned_pairs=counter[left,right,'coassigned'],
                unresolved_pairs=counter[left,right,'unresolved']))
        metrics.append(dict(pair=rep,comparison=half+'_to_full',candidate_assigned=5,one_to_one_matched_assigned=4,
                            fixed_full_DNA_cells=6,anchor_ploidy_scale=2.8,candidate_ploidy_scale=2.7+rep*.1))
    return pd.DataFrame(metrics),pd.DataFrame(pairs)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'REPORT_DESIGN.md')
    if not args.preflight:
        inputs.update(audit=args.task/'audit_v1/run_manifest.json',contrasts=args.task/'contrasts_v1/run_manifest.json',
                      metrics=args.task/'audit_v1/comparison_metrics.tsv',pairs=args.task/'contrasts_v1/cell_pair_counts.tsv',
                      helper=Path(__file__).with_name('andor_dna_stability_audit.py'))
    run=start_run(args.out,inputs,{},dict(scope='PDF display of existing descriptive diagnostics; no new scientific test'),mode='andor_dna_report')
    try:
        if args.preflight:
            metrics,pairs=fixture()
        else:
            spec=importlib.util.spec_from_file_location('verified_report',inputs['helper'])
            module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
            module.verified_run(args.task/'audit_v1');module.verified_run(args.task/'contrasts_v1')
            metrics=pd.read_csv(inputs['metrics'],sep='\t');metrics=metrics.loc[metrics.comparison.isin(['A_to_full','B_to_full'])]
            pairs=pd.read_csv(inputs['pairs'],sep='\t')
            assert len(metrics)==40 and set(metrics.pair)==set(range(20))
        render(metrics,pairs,args.out)
        outputs=sorted(args.out.glob('*.pdf'));assert len(outputs)==5
        for path in outputs:
            assert path.read_bytes().startswith(b'%PDF-')
        metrics.to_csv(args.out/'plotted_fit_values.tsv',sep='\t',index=False)
        (args.out/'captions.txt').write_text(
            'All panels: one exposed DNA specimen, 20 independent read-split seeds and 40 dependent complementary halves; '
            'no independent biological validation or new confidence claim. The rendering preflight instead uses synthetic fixtures.\n'
            'Assignment agreement: both axes retain all full-depth assay cells; agreement is with inferred DNA groups, not truth. '
            'Gray lines connect complementary halves; circles/triangles distinguish A/B.\n'
            'Integer-CN scale: gray horizontal line is the fitted full-depth reference; this scale does not identify absolute ploidy.\n'
            'Within-group coassignment: all fixed reference cell pairs remain in the denominator. A zero may reflect separation '
            'or unresolved cells; the companion unresolved panel and source tables distinguish them.\n'
            'Between-group panels: mean coassignment and unresolved fractions use fixed reference cell pairs. '
            'Only the upper triangle is displayed. Reference group numbers follow full-depth size order; '
            'the common 0-to-1 scale preserves zeros and ones. Cell pairs are not independent specimens.\n')
        write_json(dict(synthetic_render_fixture=args.preflight,biological_specimens=0 if args.preflight else 1,
                        independent_validation=False,panels=5,pdf_files=[file_record(path) for path in outputs],
                        scope='Conditional DNA measurement diagnostics, not biological clone validation',
                        visual_review='pending',new_support_cutoff=None),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
