"""Research robust intercept for the existing continuous RNA boundary model."""
import argparse
import importlib.util
from pathlib import Path
import statistics

import numpy as np
from haplodose.bulk import BULK_PROFILE
from haplodose.io import write_json


def load(name):
    spec=importlib.util.spec_from_file_location('centered_gaussian_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def checks():
    gaussian=load('reference_gaussian_regions');old=gaussian.checks()
    values=np.repeat(np.array([-1.,-.5,0.,.5,1.,2.]),2)
    expected=np.full(len(values),3.5);variance=np.ones(len(values));bins=np.repeat(np.arange(6),2)
    options=(expected,variance,bins,np.full(6,.125),np.r_[True,np.zeros(5,bool)],BULK_PROFILE)
    fitted=gaussian.fit(values,*options,centre_mode='bin-median')
    shifted=gaussian.fit(values+8.,*options,centre_mode='bin-median')
    assert fitted['centre']==statistics.median(map(float,fitted['bin_values']))
    assert shifted['centre']-fitted['centre']==8.
    for key in ('path','posterior','sampling_variance','total_variance'):np.testing.assert_array_equal(fitted[key],shifted[key])
    assert fitted['grid_scores']==shifted['grid_scores'] and fitted['beta']==shifted['beta']
    first=gaussian.fit(values,*options);second=gaussian.fit(values+8.,*options)
    assert not np.array_equal(first['posterior'],second['posterior'])
    return dict(previous_gaussian_checks=old,independent_scalar_intercept=True,
        shifted_path_posterior_and_grid_scores_exact=True,fixed_zero_counterexample=True,
        intercept_is_robust_location_not_joint_MLE=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:
        load('reference_weighted_regions').panel(argparse.Namespace(out=args.out,task=args.task,
            reference_refit=0,comparison_root=args.task.parent/'reference_gaussian_regions_20261109',
            weighting='reference',normal_centering=True,reference_gaussian_regions=True,gaussian_centre='bin-median'))


if __name__=='__main__':main()
