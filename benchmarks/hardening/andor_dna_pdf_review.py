"""Prepare temporary PDF inspection images; canonical figure artifacts remain PDFs."""
import argparse
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    tools=Path('/gpfs/commons/home/jameslee/micromamba/bin')
    helper=Path(__file__).with_name('andor_dna_stability_audit.py')
    inputs=dict(script=Path(__file__),helper=helper,report=args.task/'report_v1/run_manifest.json')
    for name in ('pdftoppm','pdftotext','pdfinfo','pdffonts'):
        inputs[name]=tools/name
    run=start_run(args.out,inputs,{},dict(scope='PDF inspection preparation only, not visual approval or scientific validation'),
                  mode='andor_dna_pdf_review')
    try:
        spec=importlib.util.spec_from_file_location('verified_pdf',helper)
        verified=importlib.util.module_from_spec(spec);spec.loader.exec_module(verified)
        verified.verified_run(args.task/'report_v1')
        environment=dict(os.environ,LD_LIBRARY_PATH=str(tools.parent/'lib'))
        version=subprocess.run([str(tools/'pdftoppm'),'-v'],capture_output=True,text=True,check=True,env=environment)
        (args.out/'poppler_version.txt').write_text(version.stdout+version.stderr)
        temporary=Path(tempfile.mkdtemp(prefix='haplodose_pdf_review_',dir='/gpfs/commons/home/jameslee/.cache'))
        records=[]
        figures=sorted((args.task/'report_v1').glob('*.pdf'));assert len(figures)==5
        for pdf in figures:
            for command in ('pdfinfo','pdffonts'):
                result=subprocess.run([str(tools/command),str(pdf)],capture_output=True,text=True,check=True,env=environment)
                (args.out/(pdf.stem+'.'+command+'.txt')).write_text(result.stdout+result.stderr)
            subprocess.run([str(tools/'pdftotext'),'-layout',str(pdf),str(args.out/(pdf.stem+'.text.txt'))],
                           check=True,capture_output=True,env=environment)
            prefix=temporary/pdf.stem
            result=subprocess.run([str(tools/'pdftoppm'),'-f','1','-singlefile','-scale-to','1800','-png',str(pdf),str(prefix)],
                                  capture_output=True,text=True,check=True,env=environment)
            assert not result.stderr.strip(),result.stderr
            image=prefix.with_suffix('.png');assert image.is_file()
            records.append(dict(PDF=file_record(pdf),temporary_inspection_image=file_record(image)))
        write_json(dict(figures=records,canonical_format='PDF',temporary_images_not_figure_deliverables=True,
                        visual_review='pending; inspect the actual images before recording a decision',
                        scientific_validation=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
