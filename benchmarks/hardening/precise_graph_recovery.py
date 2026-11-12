"""Five unchanged recovery assertions under a numerical graph correction."""
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
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mode',choices=['default','pooled_exact'],required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent
    graph=load('precision',root/'precise_proposal_graph.py');pool=load('pool',root/'pooled_initialization.py')
    script='membership_stage_panel.py' if a.mode=='default' else 'shared_membership_panel.py'
    panel=load('panel',root/script)
    inputs={name:root/name for name in ('precise_graph_recovery.py','precise_proposal_graph.py','pooled_initialization.py',
        'membership_stage_panel.py','shared_membership_panel.py','membership_weighted_stages.py','shared_membership_stages.py',
        'shared_membership_hmm.py','shared_membership_probe.py','recovery_baseline.py')}
    inputs.update(assertions=a.assertion_source,panel_manifest=a.panel_root/'manifest.json')
    rec=start_run(a.out,inputs,dict(mode=a.mode),dict(scope='Numerical precision intervention, all original recovery criteria retained'),
                  mode='precise_graph_recovery')
    trace=[];previous=sys.argv
    try:
        (a.out/'checks.json').write_text(json.dumps(graph.checks(),indent=2))
        sys.argv=[str(root/script),'--out',str(a.out/'panel'),'--panel-root',str(a.panel_root),
                  '--assertion-source',str(a.assertion_source),'--stage','hard' if a.mode=='default' else 'exact_matched']
        with graph.precise_graph(),(nullcontext() if a.mode=='default' else pool.pooled_initialization(trace)):
            panel.main()
        (a.out/'initialization_trace.json').write_text(json.dumps(trace,indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise
    finally:sys.argv=previous


if __name__=='__main__':main()
