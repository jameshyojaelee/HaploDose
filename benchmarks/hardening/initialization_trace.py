"""Read-only stage instrumentation of six locked complete-workflow replays.

No criterion or algorithm is changed. Final arrays must equal their frozen
counterparts exactly before stage traces may be interpreted.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import json
from pathlib import Path
import sys

import numpy as np
from haplodose.clones import CloneParams, SpecimenModel
from haplodose.provenance import start_run, finish_run, file_record


def load(name, path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def groups(R):
    return dict(columns=R.shape[1], map_sizes=dict(Counter(R.argmax(axis=1).tolist())),
                probability_mass=R.sum(axis=0, dtype=np.float64).tolist())


@contextmanager
def observe_phases(events, models):
    original_em, original_split = SpecimenModel.em, SpecimenModel.try_splits

    def em(self, doses, pi, *, stop_on_map=True):
        caller = sys._getframe(1)
        loc = caller.f_locals
        entry = dict(phase='em', input_groups=len(doses), cycle=loc.get('cycle'),
                     log_start=len(self.log))
        if not models:
            # fit_specimen's initialization variables, before any EM or split.
            if caller.f_code.co_name != 'fit_specimen':
                raise RuntimeError('Unexpected first EM caller')
            entry.update(initial_leiden_sizes=dict(Counter(loc['g'].tolist())),
                         initial_agglomerated_members=loc['members'],
                         initial_tumour_indices=loc['tum'])
            models.append(self)
        result = original_em(self, doses, pi, stop_on_map=stop_on_map)
        entry.update(output_groups=len(result[0]), assignments=groups(result[2]),
                     log_end=len(self.log))
        events.append(entry)
        return result

    def split(self, R, doses, pi):
        entry = dict(phase='try_splits', input_groups=len(doses),
                     assignments=groups(R), log_start=len(self.log))
        result = original_split(self, R, doses, pi)
        entry.update(output_groups=len(result[0]), accepted=bool(result[2]), log_end=len(self.log))
        events.append(entry)
        return result

    # Do not wrap merge_unsupported/agglomerate: stage hooks identify their
    # caller code objects. Read merge actions from the original model log.
    SpecimenModel.em, SpecimenModel.try_splits = em, split
    try:
        yield
    finally:
        SpecimenModel.em, SpecimenModel.try_splits = original_em, original_split


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--rank-source', type=Path, required=True)
    p.add_argument('--frozen-shards', type=Path, required=True)
    a = p.parse_args()
    runner = load('frozen_exact_runner', a.source/'exact_stage_null.py')
    rank = load('frozen_rank', a.rank_source/'global_rank_gate.py')
    gen = rank.load_generator(a.rank_source/'global_rank_generator.py')
    hook = load('frozen_exact_hooks', a.source/'shared_membership_stages.py')
    design = json.loads((a.rank_source/'prespecified_design.json').read_text())
    rank.validate_design(design)
    params = CloneParams(n_programmes=0, spatial_neighbors=0, seed=design['master_seed'], split_mode='thin')
    a.out.mkdir(exist_ok=False)
    inputs = {str(path): path for directory in (a.source, a.rank_source)
              for path in directory.glob('*.py')}
    inputs['design'] = a.rank_source/'prespecified_design.json'
    inputs['instrumentation'] = Path(__file__)
    for stage in runner.STAGES:
        for role in ('eval_nb', 'eval_state', 'eval_two'):
            inputs[f'{stage}_{role}_records'] = a.frozen_shards/stage/role/'records.jsonl'
    rec = start_run(a.out, inputs, dict(replicate=0, roles=['eval_nb','eval_state','eval_two']),
                    dict(scope='Mechanism trace; not a new performance estimate'), mode='initialization_trace')
    rows = []
    try:
        for stage in runner.STAGES:
            for role in ('eval_nb','eval_state','eval_two'):
                condition, role_id, _ = runner.ROLES[role]
                directory = a.out/f'{stage}_{role}'
                directory.mkdir()
                phases, models = [], []
                with observe_phases(phases, models):
                    row = runner.run_one(rank,gen,hook,design,params,condition,role_id,0,directory,stage)
                previous = json.loads((a.frozen_shards/stage/role/'records.jsonl').read_text().splitlines()[0])
                # JSON object keys are strings. Compare the same serialization
                # representation, while retaining bitwise array checks below.
                row = json.loads(json.dumps(row))
                if set(row) != set(previous):
                    raise AssertionError('Record schema changed')
                for key in previous:
                    if key != 'assignments' and row[key] != previous[key]:
                        raise AssertionError(f'Instrumentation changed {stage}/{role}/{key}')
                with np.load(row['assignments']['path']) as new, np.load(previous['assignments']['path']) as old:
                    if set(new.files) != set(old.files):
                        raise AssertionError('Output schema changed')
                    for name in old.files:
                        if not np.array_equal(new[name], old[name]):
                            raise AssertionError(f'Instrumentation changed final array {name}')
                model = models[0]
                result = dict(stage=stage, role=role, condition=condition, seed_key=row['seed_key'],
                              frozen_identity_passed=True, frozen_assignment=previous['assignments'],
                              replay_assignment=row['assignments'], phases=phases, model_log=model.log,
                              raw_candidate_count=row['raw_candidate_count'], score=row['score'],
                              stage_trace=row['stage_trace'], input_sha256=row['input_sha256'])
                rows.append(result)
                print(json.dumps(dict(stage=stage,role=role,identity_passed=True,
                                      groups=row['raw_candidate_count'],phases=phases)),flush=True)
        (a.out/'traces.json').write_text(json.dumps(rows, indent=2, default=lambda x:x.item()))
        finish_run(a.out, rec, backend='cpu')
    except Exception as error:
        finish_run(a.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
