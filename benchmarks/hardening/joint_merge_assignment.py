"""Research-only full-count memberships at merge in an explicitly joint fit.

This replaces neither the terminal merge rule nor the default thinning method.
Assignments and scoring share counts, so complete-workflow calibration is needed.
"""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys
import weakref
import numpy as np
from haplodose import clones,pipeline
from haplodose.provenance import start_run,finish_run


@contextmanager
def full_joint_merge(trace):
    posterior=clones.posteriors;write=pipeline.write_json;seen=weakref.WeakKeyDictionary();missing=object()
    def assign(*args,**kwargs):
        caller=sys._getframe(1)
        if caller.f_code is not clones.SpecimenModel.merge_unsupported.__code__:
            return posterior(*args,**kwargs)
        model=caller.f_locals['self']
        if model.alleles is not None or model.graph is not None or not model.thin:
            raise ValueError('Only expression-only nonspatial joint-data research mode')
        full=np.asarray(caller.f_locals['R']);proposal=posterior(*args,**kwargs)
        if full.shape!=proposal.shape or not np.isfinite(full).all():raise ValueError('Membership alignment')
        if model not in seen:seen[model]=getattr(model,'_research_merge_assignment_source',missing)
        model._research_merge_assignment_source='joint_full'
        trace.append(dict(specimen_rows=model.rows.rows.tolist(),source='joint_full',
            full_responsibility=full.tolist(),counterfactual_proposal_responsibility=proposal.tolist(),
            maximum_probability_change=float(np.max(np.abs(full-proposal)))))
        return full.copy()
    def tagged(value,path):
        if Path(path).name=='qc.json':
            value.update(merge_assignment_source='joint_full',merge_assignment_confirmation_status='not_independent',
                terminal_validation='production_merge_rule_retained')
        return write(value,path)
    clones.posteriors=assign;pipeline.write_json=tagged
    try:yield
    finally:
        clones.posteriors=posterior;pipeline.write_json=write
        for model,previous in seen.items():
            if previous is missing:del model._research_merge_assignment_source
            else:model._research_merge_assignment_source=previous


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent;script=root/'joint_unordered_panel.py'
    spec=importlib.util.spec_from_file_location('joint_panel',script);panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)
    rec=start_run(a.out,dict(script=Path(__file__),panel=script,design=a.task/'DESIGN.md'),
        dict(merge_assignment_source='joint_full'),dict(scope='Single-stage ablation with terminal correction; not independent confirmation'),
        mode='joint_merge_assignment')
    trace=[];argv=sys.argv
    try:
        sys.argv=[str(script),'--out',str(a.out/'panel'),'--task',str(a.task),'--panel-root',str(a.panel_root),
            '--assertion-source',str(a.assertion_source),'--merge','margin']
        with full_joint_merge(trace):panel.main()
        (a.out/'assignment_trace.json').write_text(json.dumps(trace,indent=2));finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        (a.out/'assignment_trace_partial.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
