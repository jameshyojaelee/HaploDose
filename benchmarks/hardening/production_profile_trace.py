"""Observe profile construction at fixed membership/nuisance without changing fits."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def capture(records, directory):
    original = clones.SpecimenModel.shared_profiles
    models = {}
    def observe(self, Yk, Lk, XLk, labels=None):
        frame = sys._getframe(1)
        assert frame.f_code is clones.SpecimenModel.em.__code__
        assert self.alleles is None and self.graph is None and self.rows.n <= 600
        if self not in models:
            name = directory/f'model_{len(models):02d}.npz'
            np.savez_compressed(name, Y=self.rows.take(np.arange(self.rows.n)), rows=self.rows.rows,
                                L=self.L, x=self.x, reset=self.reset)
            models[self] = name
        # Hold live model references to prevent Python from reusing an id.
        model_file = models[self]
        before = np.asarray(frame.f_locals['doses']).copy()
        R = frame.f_locals['R'].copy()
        fixed = {k:np.asarray(getattr(self, k)).copy() for k in ('base','offset','phi','gamma','U','W')}
        prof = original(self, Yk, Lk, XLk, labels=labels)
        assert all(np.array_equal(getattr(self,k),v) for k,v in fixed.items())
        own = prof.get('dose_own',prof['dose']).copy()
        after = prof['dose'].copy()
        own[0] = after[0] = 0.
        kernel = [self.loglik(d) for d in (before, own, after)]
        path = directory/f'profile_{len(records):04d}.npz'
        np.savez_compressed(path, R=R, before=before, own=own, after=after,
            Yk=Yk,Lk=Lk,XLk=XLk,labels=labels,centre=prof['centre'],tau=prof['tau'],
            kernel_before=kernel[0],kernel_own=kernel[1],kernel_after=kernel[2],**fixed)
        seeds = [int(p.removeprefix('seed_')) for p in Path(self.rows.Y.filename).parts if p.startswith('seed_')]
        assert len(seeds)==1
        records.append(dict(model=str(model_file), arrays=str(path), index=len(records), seed=seeds[0],
            iteration=int(frame.f_locals['iteration']), log_position=len(self.log),
            kernel_delta_own=float(np.sum(R*(kernel[1].astype(float)-kernel[0]))),
            kernel_delta_shared=float(np.sum(R*(kernel[2].astype(float)-kernel[0]))),
            kernel_delta_sharing=float(np.sum(R*(kernel[2].astype(float)-kernel[1])))))
        return prof
    clones.SpecimenModel.shared_profiles = observe
    try:
        yield
    finally:
        clones.SpecimenModel.shared_profiles = original
        models.clear()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(exist_ok=False)
    (a.out/'arrays').mkdir()
    runner=load('profile_original_recovery',Path(__file__).with_name('gene_thinning_recovery.py'))
    previous=a.task.parent/'exact_profile_merge_20261109/production_v1/panel'
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',previous=previous/'summary.json'),{},
                  dict(scope='Original5 production cases; observed fixed-membership objective, no oracle'),mode='production_profile_trace')
    records,argv=[],sys.argv
    try:
        sys.argv=[str(runner.__file__),'--out',str(a.out/'replay'),'--task',str(a.task.parent/'gene_thinning_20261109'),
                  '--mode','control','--stage','hard']
        with capture(records,a.out/'arrays'):
            runner.main()
        current=a.out/'replay/panel'
        old=json.loads((previous/'summary.json').read_text())
        new=json.loads((current/'summary.json').read_text())
        for before,after in zip(old['cases'],new['cases']):
            assert before['seed']==after['seed'] and before['criteria']==after['criteria']
            for name in ('cells.tsv.gz','clones.tsv','differences.tsv','segments.tsv','profiles.tsv.gz'):
                pd.testing.assert_frame_equal(pd.read_csv(current/f'seed_{after["seed"]}/fit'/name,sep='\t'),
                    pd.read_csv(previous/f'seed_{before["seed"]}/fit'/name,sep='\t'),check_exact=True)
            assert np.array_equal(np.load(current/f'seed_{after["seed"]}/fit/posteriors.npy'),
                                  np.load(previous/f'seed_{before["seed"]}/fit/posteriors.npy'),equal_nan=True)
        write_json(dict(all_public_outputs_identical=True,all30_criteria_identical=True,records=records,
                        production_changed=False,biological_specimens=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU observed installed NB kernel; production output parity')
    except BaseException as error:
        write_json(records,a.out/'partial_records.json')
        finish_run(a.out,rec,error=error)
        raise
    finally:
        sys.argv=argv


if __name__=='__main__':
    main()
