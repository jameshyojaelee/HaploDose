"""Audit saved DNA boundary reads without issuing further archive queries."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    source=args.task/'probe_v2'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',manifest=source/'run_manifest.json',
                queries=source/'queries.json',report=source/'report.json')
    run=start_run(args.out,inputs,{},dict(scope='Saved primary-record coordinate audit only'),mode='andor_dna_coordinate_audit')
    try:
        manifest=json.loads(inputs['manifest'].read_text());assert manifest['status']=='complete'
        for group in ('input','output'):
            for record in manifest[group]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        queries=json.loads(inputs['queries'].read_text());moved=[];all_counts=Counter();seen=set()
        for query in queries:
            assert query['exit_code']==0
            for name in ('raw','filtered','stderr'):
                assert file_record(query[name]['path'])['sha256']==query[name]['sha256']
            raw=Path(query['raw']['path']).read_text();expected=[]
            for line in raw.splitlines():
                if line.startswith('@'):
                    continue
                fields=line.split('\t');flag=int(fields[1])
                if flag&0x42==0x42 and flag&0xF04==0 and int(fields[4])>=30:
                    expected.append(line)
            observed=[line for line in Path(query['filtered']['path']).read_text().splitlines() if not line.startswith('@')]
            assert expected==observed
            boundary=int(query['region'].split(':')[1].split('-')[0])+1
            local=[];count=0
            for row in observed:
                fields=row.split('\t');tags={v.split(':',2)[0]:v.split(':',2)[2] for v in fields[11:]}
                if 'RG' not in tags:
                    continue
                count+=1;all_counts[query['line']]+=1
                pos=int(fields[3]);key=(query['accession'],fields[0],fields[1],fields[2],fields[3])
                assert key not in seen;seen.add(key)
                if pos==boundary:
                    # Boundary b is the last one-based base of [b-500000,b).
                    correct=next(j for j in range(20) if j*500000<=pos-1<(j+1)*500000)
                    wrong=next(j for j in range(21) if j*500000<=pos<(j+1)*500000)
                    assert wrong==correct+1
                    record=dict(read=fields[0],rg=tags['RG'],chromosome=fields[2],SAM_POS=pos,
                                zero_based_position=pos-1,correct_bin=correct)
                    local.append(record)
                    moved.append(dict(line=query['line'],region=query['region'],**record))
                    assert '-' not in tags['RG'] and not {'CB','CR'}&set(tags)
            assert count==query['filtered_fragments']
            assert local==query['observed_boundary_starts']
        report=json.loads(inputs['report'].read_text())
        assert len(queries)==report['successful_queries']==report['queries']==60
        assert moved==report['observed_boundary_starts'] and len(moved)==10
        for check in report['fixtures']:
            assert check['corrected_boundary_assignments']==3 and check['filtered_fragments']==8
        write_json(dict(queries_checked=len(queries),eligible_fragments=dict(all_counts),
                        boundary_starts_by_line=dict(Counter(row['line'] for row in moved)),
                        observed_boundary_starts=len(moved),all_boundary_reads_sample_index_without_cell_barcode=True,
                        biological_clone_impact='not established',historical_truth_unchanged=True),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
