"""Research initialization ablation; no production or calibration default.

Pool initial tumour groups exactly as the existing selected-null initializer,
without enabling its different split filter. Subsequent search is unchanged.
"""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys
import weakref

import numpy as np
from haplodose.clones import CloneParams, SpecimenModel
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def pooled_initialization(trace):
    original = SpecimenModel.em
    seen = weakref.WeakSet()

    def em(self, doses, pi, *, stop_on_map=True):
        if self not in seen:
            seen.add(self)
            loc = sys._getframe(1).f_locals
            if sys._getframe(1).f_code.co_name != 'fit_specimen':
                raise RuntimeError('Unexpected initial EM context')
            if self.alleles is not None or self.graph is not None or self.params.split_null_draws:
                raise ValueError('Only expression-only nonspatial unfilted initialization supported')
            entry = dict(input_groups=len(doses), changed=False, rows=self.rows.rows.tolist())
            if len(doses) > 2:
                tum = loc['tum']
                if len(tum) != len(doses)-1:
                    raise RuntimeError('Initial tumour correspondence changed')
                Y, L, XL = loc['Yk'][tum].sum(axis=0), loc['Lk'][tum].sum(), loc['XLk'][tum].sum(axis=0)
                dose = self.profiles(Y, L, XL)['dose'][0]
                pi = np.array([pi[0], pi[1:].sum()])
                doses = [doses[0], dose]
                self.log.append(dict(action='experimental_pool_initial', groups=len(tum),
                                     reason='initialize one RNA tumour profile before unchanged split search'))
                entry.update(changed=True, pooled_counts=Y.tolist(), pooled_depth=float(L),
                             pooled_depth_moments=XL.tolist(), mixture_weights=pi.tolist())
            entry['output_groups'] = len(doses)
            trace.append(entry)
        return original(self, doses, pi, stop_on_map=stop_on_map)

    SpecimenModel.em = em
    try:
        yield
    finally:
        SpecimenModel.em = original


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--mode', choices=['mechanism', 'recovery'], required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--rank-source', type=Path)
    p.add_argument('--frozen-shards', type=Path)
    p.add_argument('--panel-root', type=Path)
    p.add_argument('--assertion-source', type=Path)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    inputs = {'intervention': Path(__file__)}
    for directory in (a.source, root):
        for name in ('shared_membership_panel.py','shared_membership_stages.py','shared_membership_hmm.py',
                     'membership_weighted_stages.py','exact_stage_null.py','recovery_baseline.py'):
            inputs[str(directory/name)] = directory/name
    if a.assertion_source: inputs['assertions'] = a.assertion_source
    if a.panel_root: inputs['panel_manifest'] = a.panel_root/'manifest.json'
    if a.rank_source:
        for name in ('global_rank_gate.py','global_rank_generator.py','prespecified_design.json'):
            inputs[name] = a.rank_source/name
    if a.frozen_shards:
        for role in ('eval_nb','eval_state','eval_two'):
            inputs['frozen_'+role] = a.frozen_shards/'exact_matched'/role/'records.jsonl'
    rec = start_run(a.out, inputs, dict(mode=a.mode,stage='exact_matched',criteria='unchanged'),
                    dict(scope='Initialization-only development ablation; no error or biological guarantee'),
                    mode='pooled_initialization')
    trace = []
    try:
        with pooled_initialization(trace):
            if a.mode == 'recovery':
                panel = load('locked_panel', a.source/'shared_membership_panel.py')
                old_argv = sys.argv
                sys.argv = [str(a.source/'shared_membership_panel.py'),'--out',str(a.out/'panel'),
                            '--panel-root',str(a.panel_root),'--assertion-source',str(a.assertion_source),
                            '--stage','exact_matched']
                try:
                    panel.main()
                finally:
                    sys.argv = old_argv
            else:
                runner = load('locked_runner', a.source/'exact_stage_null.py')
                rank = load('locked_rank', a.rank_source/'global_rank_gate.py')
                gen = rank.load_generator(a.rank_source/'global_rank_generator.py')
                hook = load('locked_hook', a.source/'shared_membership_stages.py')
                design = json.loads((a.rank_source/'prespecified_design.json').read_text())
                rank.validate_design(design)
                params = CloneParams(n_programmes=0,spatial_neighbors=0,seed=design['master_seed'],split_mode='thin')
                records = []
                for role in ('eval_nb','eval_state','eval_two'):
                    condition, role_id, _ = runner.ROLES[role]
                    directory = a.out/role
                    directory.mkdir()
                    row = runner.run_one(rank,gen,hook,design,params,condition,role_id,0,directory,'exact_matched')
                    old = json.loads((a.frozen_shards/'exact_matched'/role/'records.jsonl').read_text().splitlines()[0])
                    if row['input_sha256'] != old['input_sha256']:
                        raise AssertionError('Input changed')
                    if role == 'eval_nb':
                        with np.load(row['assignments']['path']) as new, np.load(old['assignments']['path']) as prior:
                            for name in prior.files:
                                np.testing.assert_array_equal(new[name],prior[name])
                        if trace[-1]['changed']: raise AssertionError('One-group control was modified')
                    row['previous_raw_candidate_count'] = old['raw_candidate_count']
                    records.append(row)
                    print(json.dumps(dict(role=role,groups=row['raw_candidate_count'],
                                          previous_groups=old['raw_candidate_count'],
                                          ARI=row['ARI_against_simulated_CNA_labels'])),flush=True)
                (a.out/'records.json').write_text(json.dumps(records,indent=2))
        (a.out/'initialization_trace.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        (a.out/'initialization_trace_partial.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,error=error)
        raise


if __name__ == '__main__':
    main()
