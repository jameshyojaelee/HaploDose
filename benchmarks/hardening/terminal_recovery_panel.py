"""Unchanged five-case recovery tracks under the terminal-validation wheel."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import sys
from haplodose.provenance import start_run,finish_run


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--mode',choices=['original','pooled_exact','joint'],required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent
    script={'original':'membership_stage_panel.py','pooled_exact':'shared_membership_panel.py','joint':'joint_unordered_panel.py'}[a.mode]
    panel=load('terminal_panel',root/script);pool=load('terminal_pool',root/'pooled_initialization.py')
    inputs={name:root/name for name in ('terminal_recovery_panel.py','membership_stage_panel.py','shared_membership_panel.py',
        'joint_unordered_panel.py','joint_unordered_stages.py','shared_membership_stages.py','membership_weighted_stages.py',
        'shared_membership_hmm.py','shared_membership_probe.py','label_invariant_membership.py','label_invariance_probe.py',
        'pooled_initialization.py','recovery_baseline.py')}
    inputs.update(design=a.task/'DESIGN.md',panel_manifest=a.panel_root/'manifest.json',assertions=a.assertion_source)
    rec=start_run(a.out,inputs,dict(mode=a.mode),dict(scope='Mandatory stopping correction; original thresholds/failures retained'),
        mode='terminal_recovery')
    trace=[];argv=sys.argv
    try:
        sys.argv=[str(root/script),'--out',str(a.out/'panel'),'--panel-root',str(a.panel_root),'--assertion-source',str(a.assertion_source)]
        if a.mode=='joint':sys.argv+=['--task',str(a.task),'--merge','margin']
        else:sys.argv+=['--stage','hard' if a.mode=='original' else 'exact_matched']
        with (pool.pooled_initialization(trace) if a.mode=='pooled_exact' else nullcontext()):panel.main()
        (a.out/'initialization_trace.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise
    finally:sys.argv=argv


if __name__=='__main__':main()
