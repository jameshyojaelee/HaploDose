"""Recheck the independently optimized saved thinning proposal on a candidate wheel."""
import argparse
import json
from pathlib import Path

import numpy as np

from haplodose.clones import _fit_event_mixture
from haplodose.io import new_output_dir, write_json
from haplodose.provenance import start_run, finish_run


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    out = new_output_dir(args.out)
    inputs = dict(runner=Path(__file__), reference=args.reference / 'summary.json',
                  arrays=args.reference / 'independent_region.npz')
    rec = start_run(out, inputs, {}, {}, mode='saved_thinning_counterexample')
    try:
        ref = json.loads(inputs['reference'].read_text())
        with np.load(inputs['arrays']) as z:
            fit = _fit_event_mixture(z['O'], z['E'], z['Rnb'], .19881382302316983)
            labels = fit['responsibility'] > .5
            agreement = max(float(np.mean(labels == z['labels'])), float(np.mean(labels != z['labels'])))
        write_json(dict(fit=fit, original_gain=ref['original_gain'], independent_gain=ref['best']['gain'],
                        gain_difference=fit['gain'] - ref['best']['gain'],
                        label_agreement=agreement, scope='Identical saved proposal counts and expectations; no threshold change'),
                   out / 'summary.json')
        finish_run(out, rec, backend='cpu')
    except Exception as error:
        finish_run(out, rec, backend='cpu', error=error)
        raise


if __name__ == '__main__':
    main()
