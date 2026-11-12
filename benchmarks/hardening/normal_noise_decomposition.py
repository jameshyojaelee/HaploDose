"""Research subtraction of approximate normal sampling variance."""
import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from haplodose.io import write_json


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('normal_noise_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def component_variance(observed,sampling):
    v,s=np.broadcast_arrays(np.asarray(observed,float),np.asarray(sampling,float))
    if not np.isfinite(v).all() or not np.isfinite(s).all() or np.any(v<0) or np.any(s<0):
        raise ValueError('Finite nonnegative variance components required')
    return np.maximum(v-s,0.)


def decompose(normal_effects,normal_expected,gene_bins,training):
    gaussian=load('reference_gaussian_regions');weighted=load('reference_weighted_regions')
    effects,expected=np.asarray(normal_effects,float),np.asarray(normal_expected,float)
    rows=[]
    for i in range(len(effects)):
        other=np.delete(effects,i,axis=0)
        adjusted,centre,_=weighted.normal_centered_effect(effects[i],other,gene_bins)
        spread=(np.median(abs(other-centre),axis=0)/ndtri(.75))**2
        variance=weighted.gene_variance(spread,expected[i],'reference')
        values,sampling=gaussian.bin_observations(adjusted,expected[i],variance,gene_bins)
        np.testing.assert_array_equal(values,training['normal_bins'][i]);rows.append(sampling)
    rows=np.asarray(rows);typical=np.median(rows,axis=0)
    intrinsic=component_variance(training['reference_variance'],typical)
    return dict(normal_sampling_variance=rows,typical_normal_sampling_variance=typical,
        observed_reference_variance=training['reference_variance'],intrinsic_variance=intrinsic)


def checks():
    observed=np.array([1.,.125,.75]);sampling=np.array([.25,.25,0.])
    np.testing.assert_array_equal(component_variance(observed,sampling),[.75,0.,.75])
    np.testing.assert_array_equal(component_variance(observed,np.zeros(3)),observed)
    # Exact binary-rational component fixture: latent3/4, normal noise1/4,
    # query noise1/2. Additive normal/query treatment gives3/2 instead of5/4.
    assert component_variance(1.,.25)+.5==1.25
    invalid=0
    for v,s in [(-1.,0.),(1.,-1.),(np.nan,0.)]:
        try:component_variance(v,s)
        except ValueError:invalid+=1
    assert invalid==3
    normals=np.array([[0.,1.,2.,3.],[1.,2.,4.,5.],[2.,3.,6.,7.]])
    expected=np.full_like(normals,3.5);bins=np.array([0,0,1,1]);gaussian=load('reference_gaussian_regions')
    training=gaussian.train(normals,expected,bins);result=decompose(normals,expected,bins,training)
    permutation=np.array([3,0,2,1]);alternate=decompose(normals[:,permutation],expected[:,permutation],bins[permutation],
        gaussian.train(normals[:,permutation],expected[:,permutation],bins[permutation]))
    for key in result:np.testing.assert_array_equal(result[key],alternate[key])
    assert np.all(result['intrinsic_variance']>=0)
    assert np.all(result['intrinsic_variance']<=result['observed_reference_variance'])
    return dict(exact_binary_rational_variance_identity=True,nonnegative_boundary=True,no_noise_limit_exact=True,
        all_normal_summaries_exact=True,gene_permutation_exact=True,invalid_components_rejected=invalid,
        original_gaussian_checks=gaussian.checks(),sampling_variance_calibrated=False)


def effect_summary(args):
    root=args.task.parents[1];new=args.task/'bulk_v1';old=root/'dev/reference_gaussian_regions_20261109/bulk_v1'
    with np.load(new/'normal_noise_decomposition.npz') as z:
        observed=z['observed_reference_variance'];intrinsic=z['intrinsic_variance']
        normal=dict(bins=len(observed),changed_bins=int(np.sum(observed!=intrinsic)),
            zero_intrinsic_bins=int(np.sum(intrinsic==0)),maximum_removed_variance=float(np.max(observed-intrinsic)),
            maximum_normal_sampling_variance=float(np.max(z['typical_normal_sampling_variance'])))
    contexts=[]
    for path in sorted(new.glob('Gaussian_context_*.npz')):
        with np.load(path) as a,np.load(old/path.name) as b:
            contexts.append(dict(file=path.name,path_exact=bool(np.array_equal(a['path'],b['path'])),
                posterior_exact=bool(np.array_equal(a['posterior'],b['posterior'])),
                changed_variance_bins=int(np.sum(a['total_variance']!=b['total_variance'])),
                maximum_variance_change=float(np.max(abs(a['total_variance']-b['total_variance'])))))
    assert len(contexts)==48
    exports=[]
    with np.load(new/'RNA_predictions.npz') as a,np.load(old/'RNA_predictions.npz') as b:
        for name in a.files:
            previous=name.replace('decomposed_gaussian_','gaussian_region_')
            exports.append(dict(method=name,preceding=previous,exact=bool(np.array_equal(a[name],b[previous],equal_nan=True))))
    first=json.loads((new/'Gaussian_decisions.json').read_text());second=json.loads((old/'Gaussian_decisions.json').read_text())
    assert [r['sample'] for r in first]==[r['sample'] for r in second]
    write_json(dict(normal=normal,contexts=contexts,exports=exports,
        changed_beta=sum(a['beta']!=b['beta'] for a,b in zip(first,second)),
        changed_tau=sum(a['tau']!=b['tau'] for a,b in zip(first,second)),
        primary_gate_changed=False,larger_uncertainty_run=False),args.out/'report.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--effect-summary',action='store_true')
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    elif args.effect_summary:effect_summary(args)
    else:
        load('reference_weighted_regions').panel(argparse.Namespace(out=args.out,task=args.task,
            reference_refit=0,comparison_root=args.task.parent/'centered_gaussian_regions_20261109',
            weighting='reference',normal_centering=True,reference_gaussian_regions=True,gaussian_variance='noise-subtracted'))


if __name__=='__main__':main()
