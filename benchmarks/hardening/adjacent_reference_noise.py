"""Research reference-only NB moments from disjoint adjacent-depth pairs."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('adjacent_reference_'+name,Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def estimate(counts, exposure, groups, cell_ids):
    y,L = np.asarray(counts,float),np.asarray(exposure,float)
    groups,ids = np.asarray(groups),np.asarray(cell_ids)
    if (y.ndim != 2 or L.shape != (len(y),) or groups.shape != L.shape or ids.shape != L.shape
            or not pd.Index(ids).is_unique or pd.Index(ids).hasnans or pd.Index(groups).hasnans
            or not np.isfinite(y).all() or np.any(y < 0) or np.any(y != np.floor(y))
            or not np.isfinite(L).all() or np.any(L <= 0)):
        raise ValueError('Identified reference counts, positive totals and declared groups required')
    if not np.array_equal(y.sum(axis=1),L):
        raise ValueError('Exposure must equal the selected-feature library total')
    moments,details,total_pairs = [],[],0
    for group in sorted(set(groups)):
        use = np.flatnonzero(groups == group)
        use = use[np.lexsort((ids[use],L[use]))]
        if len(use) < 2:
            raise ValueError('Every fit-reference group requires a pair; no fallback')
        pairs = use[:2*(len(use)//2)].reshape(-1,2)
        a,b = pairs.T
        first = (y[a]*(y[a]-1))*(L[b]/L[a])[:,None]
        second = (y[b]*(y[b]-1))*(L[a]/L[b])[:,None]
        cross = 2*y[a]*y[b]
        A = math.fsum([*map(float,first.ravel()),*map(float,second.ravel()),
                       *map(float,-cross.ravel())])
        B = math.fsum(map(float,cross.ravel()))
        magnitude = math.fsum([*map(float,first.ravel()),*map(float,second.ravel()),
                               *map(float,cross.ravel())])
        moments.append((A,B,magnitude))
        total_pairs += len(pairs)
        details.append(dict(group=str(group),reference_cells=len(use),pairs=pairs.tolist(),
            paired_cell_ids=ids[pairs].tolist(),unpaired_cell_ids=ids[use[len(pairs)*2:]].tolist(),
            numerator=A,denominator=B,magnitude=magnitude))
    numerator = math.fsum(m[0] for m in moments)
    denominator = math.fsum(m[1] for m in moments)
    magnitude = math.fsum(m[2] for m in moments)
    operations = 32*total_pairs*y.shape[1]+128
    eps = np.finfo(float).eps
    if operations*eps >= 1:
        raise ArithmeticError('Moment rounding budget cannot be represented')
    gamma = operations*eps/(1-operations*eps)
    error = gamma*max(1.,magnitude)
    if denominator <= error or numerator <= error:
        raise ValueError('Positive cross-cell signal and dispersion not numerically identified; no floor')
    phi = numerator/denominator
    phi_error = error/denominator+abs(numerator)*error/(denominator*(denominator-error))
    if phi <= phi_error:
        raise ValueError('Positive dispersion unresolved beyond arithmetic error')
    return dict(estimated_phi=phi,estimated_size=1/phi,numerator=numerator,denominator=denominator,
        numerical_moment_error=error,numerical_phi_error=phi_error,
        numerical_size_error=phi_error/(phi*(phi-phi_error)),groups=details,pairs=total_pairs,
        interpretation='Adjacent-depth reference RNA residual dispersion; not identified pure technical noise',
        calibration_status='unvalidated',uncertainty_units='arithmetic_error_only_not_sampling_SE')


def estimator(records, depth_strata=5):
    if depth_strata != 5:
        raise ValueError('Use with the unchanged five-quintile comparator workflow')
    cache = {}
    def size(info):
        frame = sys._getframe(1)
        while frame is not None and frame.f_code is not pipeline._run_single.__code__:
            frame = frame.f_back
        if frame is None:
            raise RuntimeError('Missing inference input correspondence')
        loc = frame.f_locals
        used,meta = loc['used'],loc['meta']
        key = (str(loc['counts']),tuple(used.column))
        if key in cache:
            return cache[key]
        rows = np.flatnonzero(loc['fit'])
        assert len(rows) >= 20 and meta.iloc[rows].role.eq('reference').all()
        assert meta.iloc[rows].ref_split.eq('fit').all()
        X,gbin = pipeline.gene_loader(loc['counts'],used,rows)(np.arange(len(rows)))
        np.testing.assert_array_equal(gbin,info['gbin'])
        L = np.asarray(loc['Lbin'][rows],float)
        groups = (meta.dataset+'|'+meta.ref_group).to_numpy()[rows]
        result = estimate(X.toarray(),L,groups,meta.cell_id.to_numpy()[rows])
        records.append(dict(counts=str(loc['counts']),training_rows=rows.tolist(),
            training_cell_ids=meta.cell_id.to_numpy()[rows].tolist(),genes=X.shape[1],
            ordered_genes=used.gene.tolist(),query_rows_used=0,holdout_reference_rows_used=0,
            true_noise_used=False,method='disjoint_adjacent_depth_pairs',**result))
        cache[key] = result['estimated_size']
        return cache[key]
    return size


def checks():
    y = np.array([[1,9],[8,2],[3,7],[6,4]])
    L,groups,ids = np.full(4,10.),np.array(['a']*4),np.array(['a','b','c','d'])
    fitted = estimate(y,L,groups,ids)
    phi = float(Fraction(19,36))
    assert abs(fitted['estimated_phi']-phi) <= fitted['numerical_phi_error']
    assert fitted['numerator'] == 76 and fitted['denominator'] == 144
    order = np.array([3,1,0,2])
    perm = estimate(y[order],L[order],groups[order],ids[order])
    assert perm['estimated_phi'] == fitted['estimated_phi']
    assert perm['groups'][0]['paired_cell_ids'] == fitted['groups'][0]['paired_cell_ids']
    genes = estimate(y[:,::-1],L,groups,ids)
    assert genes['estimated_phi'] == fitted['estimated_phi']
    odd = estimate(np.vstack([y,[4,6]]),np.full(5,10.),np.array(['a']*5),np.array(['a','b','c','d','e']))
    assert odd['estimated_phi'] == fitted['estimated_phi']
    assert odd['groups'][0]['unpaired_cell_ids'] == ['e']
    expectations = []
    for first,second in ((Fraction(1,4),Fraction(1,4)),(Fraction(1,4),Fraction(3,8))):
        li,lj,dispersion = Fraction(8),Fraction(12),Fraction(1,10)
        Fi,Fj = (1+dispersion)*(li*first)**2,(1+dispersion)*(lj*second)**2
        cross = li*lj*first*second
        A = lj/li*Fi+li/lj*Fj-2*cross
        B = 2*cross
        bias = (1+dispersion)*li*lj*(first-second)**2
        assert A-dispersion*B == bias
        expectations.append(dict(numerator=str(A),denominator=str(B),nonconstant_rate_bias=str(bias)))
    rejected = 0
    for counts,exposure,group,cell in ((y,L,['a','b','c','d'],ids),
        (np.ones((4,2)),np.full(4,2.),groups,ids),
        (y,L+1,groups,ids),(y,L,groups,['a']*4),
        (np.array([[0,10],[10,0]]),[10,10],['a','a'],['a','b'])):
        try:estimate(counts,exposure,group,cell)
        except ValueError:rejected += 1
        else:raise AssertionError('Invalid or unidentified reference must reject')
    return dict(exact_sample_phi='19/36',sample_phi_error=abs(fitted['estimated_phi']-phi),
        numerical_phi_budget=fitted['numerical_phi_error'],exact_expectations=expectations,
        canonical_depth_ties=True,row_permutation=True,gene_permutation=True,odd_unused_cell_explicit=True,
        invalid_or_unidentified_rejected=rejected,
        biological_calibration=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--numerics',action='store_true')
    parser.add_argument('--control',action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for name in ('gene_shared_recovery','gene_shared_confirmation','reference_gene_noise','conditioned_gene_score',
            'gene_thinning_hook','depth_conditioned_genes','membership_weighted_stages','gene_veto_weights'):
        inputs[name] = Path(__file__).with_name(name+'.py')
    if not args.numerics:
        inputs['original_manifest'] = args.task.parent/'recovery_20261107/recovery_panel_v1/manifest.json'
    rec = start_run(args.out,inputs,dict(control=args.control,numerics=args.numerics),
        dict(scope='Reference-noise method and original recovery; no biological or5% claim'),mode='adjacent_reference_noise')
    try:
        if args.numerics:
            write_json(checks(),args.out/'report.json')
        else:
            runner = load('gene_shared_recovery')
            original_load,argv = runner.load,sys.argv
            def altered(name):
                module = original_load(name)
                if name == 'reference_gene_noise' and not args.control:
                    module.estimator = estimator
                return module
            runner.load = altered
            try:
                sys.argv = [str(runner.__file__),'--task',str(args.task),'--out',str(args.out/'workflow'),
                            '--mode','depth_gene_score']
                runner.main()
            finally:
                runner.load,sys.argv = original_load,argv
            report = json.loads((args.out/'workflow/report.json').read_text())
            exact = None
            if args.control:
                old = args.task.parent/'conditioned_gene_score_20261109/full_candidate_v2'
                previous = json.loads((old/'report.json').read_text())
                assert previous['recovered'] == report['recovered'] == 5
                assert [r['criteria'] for r in previous['cases']] == [r['criteria'] for r in report['cases']]
                arrays,tables,folds = 0,0,0
                for row in report['cases']:
                    before = old/'candidate/panel'/f"seed_{row['seed']}"/'fit'
                    after = args.out/'workflow/candidate/panel'/f"seed_{row['seed']}"/'fit'
                    for path in before.iterdir():
                        if path.suffix == '.npy':
                            np.testing.assert_array_equal(np.load(path),np.load(after/path.name));arrays += 1
                        elif path.name.endswith(('.tsv','.tsv.gz')):
                            pd.testing.assert_frame_equal(pd.read_csv(path,sep='\t'),pd.read_csv(after/path.name,sep='\t'),check_exact=True);tables += 1
                for path in (old/'candidate/folds').glob('*.npz'):
                    with np.load(path) as before,np.load(args.out/'workflow/candidate/folds'/path.name) as after:
                        assert before.files == after.files
                        for key in before.files:np.testing.assert_array_equal(before[key],after[key])
                    folds += 1
                assert folds == 15
                exact = dict(criteria=30,arrays=arrays,tables=tables,folds=folds,exact=True)
            write_json(dict(control=args.control,recovered=report['recovered'],cases=report['cases'],
                strict_control=exact,production_changed=False,scientific_criteria_changed=False,
                candidate_passed_all_original=bool(report['recovered']==5)),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU reference-only paired-depth moments and unchanged full gene-count recovery')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__ == '__main__':
    main()
