"""Preflight independent DNA metrics before corrected real outputs are inspected."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    reference=Path(__file__).with_name('andor_dna_independent.py')
    spec=importlib.util.spec_from_file_location('independent_dna',reference)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=args.task/'preflight_v3'
    run=start_run(args.out,dict(script=Path(__file__),reference=reference,
                  design=args.task/'AUDIT_DESIGN.md',fixture=source/'run_manifest.json',
                  expected=source/'report.json'),{},dict(scope='Independent numerical references, biological n0'),
                  mode='andor_dna_independent_preflight')
    try:
        same=module.partition_reference([0,0,1,1],[9,9,8,8]);assert same['matched']==4 and same['ARI']==1.
        crossed=module.partition_reference([0,0,1,1],[0,1,0,1]);assert crossed['matched']==2 and crossed['ARI']==-.5
        rng=np.random.default_rng([20261108,1410]);exhaustive=0
        for nrow,ncol in itertools.product(range(1,4),range(1,5)):
            for rep in range(5):
                table=rng.integers(0,5,(nrow,ncol));left=[];right=[]
                for i,j in itertools.product(range(nrow),range(ncol)):
                    left.extend([i]*int(table[i,j]));right.extend([j]*int(table[i,j]))
                result=module.partition_reference(left,right)
                if nrow<=ncol:
                    best=max(sum(int(table[i,j]) for i,j in enumerate(p)) for p in itertools.permutations(range(ncol),nrow))
                else:
                    best=max(sum(int(table[i,j]) for j,i in enumerate(p)) for p in itertools.permutations(range(nrow),ncol))
                assert result['matched']==best;exhaustive+=1
        expected=json.loads((source/'report.json').read_text())
        audits=[]
        for left,right,name in (('full','full','self_comparison'),('A','B','half_comparison')):
            audits.append(module.audit_comparison(source/left/'DNA/SNU-601',source/right/'DNA/SNU-601',expected[name]))
        corrupted=dict(expected['half_comparison'],both_assigned=48)
        try:
            module.audit_comparison(source/'A/DNA/SNU-601',source/'B/DNA/SNU-601',corrupted)
        except AssertionError:
            pass
        else:
            raise AssertionError('Corrupt assignment denominator accepted')
        # Exact weighted 30% eligibility, independent of a convenient bin count.
        for name,length,qualifies in (('exact',3,True),('below',2,False)):
            directory=args.out/name;directory.mkdir()
            write_json(dict(clones={'dA':1,'dB':1}),directory/'qc.json')
            (directory/'clone_bins.tsv').write_text('chromosome\tstart\tend\tarm\tdA\tdB\n'
                f'chr1\t0\t{length}\tchr1q\t2\t1\nchr1\t{length}\t10\tchr1q\t2\t2\n')
            (directory/'clone_differences.tsv').write_text('clone_a\tclone_b\tarm\tfrac_bins_diff\tmean_diff\n'
                +('dA\tdB\tchr1q\t0.3\t1\n' if qualifies else ''))
            assert module.arm_reference(directory)==(['chr1q'] if qualifies else [])
        write_json(dict(exhaustive_matching_tables=exhaustive,hand_partition_ARI_exact=True,
                        synthetic_artifact_comparisons=audits,corrupted_denominator_rejected=True,
                        supported_bp_30percent_boundary=True,biological_specimens=0),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
