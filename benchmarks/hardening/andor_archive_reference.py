"""Bounded archive metadata investigation after header-only alignment queries."""
import argparse
from pathlib import Path
import subprocess

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),{},
                  dict(scope='Archive reference metadata only; no truth reconstruction'),mode='andor_archive_reference')
    try:
        root=Path('/gpfs/commons/home/jameslee/Haplodose/data/external/andor2020_scdna/sra')
        records=[]
        commands=[('sam_dump_help',['sam-dump','--help'])]
        for accession in ('SRR8117283','SRR8117281','SRR8117282'):
            path=root/accession/(accession+'.sra')
            commands.append((accession+'_references',['vdb-dump',str(path),'-T','REFERENCE','-C','NAME,SEQ_ID','-R','1-3']))
            commands.append((accession+'_alignment',['vdb-dump',str(path),'-T','PRIMARY_ALIGNMENT','-C','REF_ID,REF_POS,SAM_FLAGS','-R','1-3']))
        # The first observed SNU-601 alignment has zero-based REF_POS9991.
        # A known-read control distinguishes query behavior from no coverage.
        snu=root/'SRR8117283/SRR8117283.sra'
        for label,region in [('known_point','chr1:9992-9992'),('known_window','chr1:9500-10500'),
                             ('boundary_window','chr1:499500-500500')]:
            commands.append((label,['sam-dump','--primary','--aligned-region',region,str(snu)]))
        for label,command in commands:
            result=subprocess.run(command,capture_output=True,text=True,timeout=120)
            output=args.out/(label+'.stdout');output.write_text(result.stdout)
            error=args.out/(label+'.stderr');error.write_text(result.stderr)
            records.append(dict(label=label,command=command,exit_code=result.returncode,
                                stdout=file_record(output),stderr=file_record(error)))
        write_json(records,args.out/'commands.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
