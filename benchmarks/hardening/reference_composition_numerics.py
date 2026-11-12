"""Independent reference-composition proportions and transport invariants."""
import argparse
from fractions import Fraction
import importlib.util
import math
from pathlib import Path

import numpy as np
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('reference_gene_composition.py')
    spec = importlib.util.spec_from_file_location('checked_reference_composition',path)
    method = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(method)
    run = start_run(args.out,dict(script=Path(__file__),method=path,design=args.task/'DESIGN.md'),{},
        dict(scope='Numerical references, not biological validation'),mode='reference_composition_numerics')
    try:
        X = np.array([[2,6,2,4],[2,0,4,0],[1,3,2,6],[3,1,6,2]],float)
        depth = X.sum(axis=1)
        groups,gbin = np.array(['A','A','B','B']),np.array([0,0,1,1])
        fitted = method.table(X,depth,groups,gbin,strata=1)
        expected = np.array([float(Fraction(19,44)),float(Fraction(25,44)),
                             float(Fraction(23,44)),float(Fraction(21,44))])
        predicted = method.predict(fitted,[10.,20.],gbin,[[3.,1.],[1.,3.]])
        operations = 8*X.size+32
        unit = np.finfo(float).eps
        budget = operations*unit/(1-operations*unit)
        np.testing.assert_allclose(predicted,np.broadcast_to(expected,predicted.shape),rtol=0.,atol=budget)
        for b in range(2):
            np.testing.assert_allclose(predicted[:,gbin == b].sum(axis=1),1.,rtol=0.,atol=budget)
        row_order,gene_order = np.array([3,0,2,1]),np.array([2,0,3,1])
        reordered = method.table(X[row_order][:,gene_order],depth[row_order],groups[row_order],gbin[gene_order],strata=1)
        other = method.predict(reordered,[10.,20.],gbin[gene_order],[[3.,1.],[1.,3.]])
        np.testing.assert_allclose(other,predicted[:,gene_order],rtol=0.,atol=budget)
        counts = np.array([[i,i*i+3] for i in range(1,11)],float)
        L = counts.sum(axis=1)
        bins = np.zeros(2,int)
        curve = method.table(counts,L,np.repeat('A',10),bins,strata=5)
        independent_fractions,independent_knots = [],[]
        for s in range(5):
            rows = [2*s,2*s+1]
            exposure = sum(int(L[c]) for c in rows)
            weight = Fraction(exposure,int(L.sum()))
            totals = [sum(int(counts[c,g]) for c in rows)+Fraction(1,2)*weight for g in range(2)]
            independent_fractions.append([float(x/sum(totals)) for x in totals])
            independent_knots.append(math.fsum(math.log(float(L[c]))*float(L[c]) for c in rows)/exposure)
        np.testing.assert_allclose(curve['fraction'][0],independent_fractions,rtol=0.,atol=budget)
        np.testing.assert_allclose(curve['knots'][0],independent_knots,rtol=0.,atol=budget*max(independent_knots))
        query = np.exp((independent_knots[1]+independent_knots[2])/2)
        interpolated = method.predict(curve,[query],bins,[[1.]])[0]
        scalar = (np.array(independent_fractions[1])+np.array(independent_fractions[2]))/2
        np.testing.assert_allclose(interpolated,scalar,rtol=0.,atol=budget)
        tails = method.predict(curve,[L.min()/10,L.max()*10],bins,[[1.]])
        np.testing.assert_allclose(tails,[independent_fractions[0],independent_fractions[-1]],rtol=0.,atol=budget)
        try:
            method.table(X,depth,groups,gbin,strata=5)
        except ValueError:
            pass
        else:
            raise AssertionError('Unsupported reference strata must fail')
        write_json(dict(all_checks_passed=True,Fraction_mixed_reference_proportions=True,
            row_gene_permutation=True,independent_depth_proportions_and_knots=True,
            scalar_interpolation_and_endpoint_rule=True,unsupported_reference_abstention=True,
            operation_rounding_allowance=budget),args.out/'report.json')
        finish_run(args.out,run,backend='CPU independent rational reference proportions and scalar interpolation')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
