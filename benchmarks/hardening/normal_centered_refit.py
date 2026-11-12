"""Reconstruct complete comparator families before normal-centered RNA profiles."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def load(name):
    spec=importlib.util.spec_from_file_location('normal_centered_refit_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--reference-refit',type=int,choices=range(9),required=True);args=parser.parse_args()
    args.out.mkdir(exist_ok=False);number=args.reference_refit
    comparison=args.task/'comparators';comparison.mkdir(exist_ok=True)
    previous=comparison/('bulk_v1' if number==0 else f'refit_{number:03d}_v1')
    saved=sys.argv
    try:
        selector=load('rna_weight_selection')
        # Keep its original design/provenance while reconstructing the same
        # failed comparator on each existing reference; no re-selection by DNA.
        sys.argv=[str(selector.__file__),'--out',str(previous),'--task',str(args.task.parent/'rna_weight_selection_20261109'),
            '--reference-refit',str(number)]
        selector.main()
        weighted=load('reference_weighted_regions')
        weighted.panel(argparse.Namespace(out=args.out,task=args.task,reference_refit=number,
            comparison_root=comparison,weighting='reference',normal_centering=True))
        if number==0:
            prior=args.task/'bulk_v1'
            with np.load(prior/'RNA_predictions.npz') as old,np.load(args.out/'RNA_predictions.npz') as new:
                assert old.files==new.files
                for key in old.files:np.testing.assert_array_equal(old[key],new[key])
            pd.testing.assert_frame_equal(pd.read_csv(prior/'patient_scores.tsv',sep='\t',float_precision='round_trip'),
                pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip'),check_exact=True)
        (args.out/'refit_control.json').write_text(json.dumps(dict(reference_refit=number,
            preceding_methods_retained=24,primary_exact=(number==0),criteria_changed=False),indent=2)+'\n')
    finally:sys.argv=saved


if __name__=='__main__':main()
