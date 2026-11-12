"""Reconstruct old/corrected DNA in one environment before conditional resampling."""
import argparse
import importlib.util
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('andor_dna_stability_common.py')
    spec=importlib.util.spec_from_file_location('dna_common',helper)
    common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
    common.verify_run(args.task/'count_audit_v1')
    module=common.dna_module()
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    old=root/'data/external/andor2020_scdna/bins_original_bam'
    new=root/'dev/andor_dna_coordinates_20261108/corrected_counts_v1/bins'
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,DNA_source=Path(module.__file__),
                  design=args.task/'DESIGN.md',count_gate=args.task/'count_audit_v1/run_manifest.json',
                  old=old/'SNU-601.tsv.gz',new=new/'SNU-601.tsv.gz'),{},
                  dict(scope='One exposed DNA specimen, identical parameters/environment, no RNA input'),mode='andor_dna_reconstruction')
    try:
        old_fit=common.fit(module,old,args.out/'old')
        new_fit=common.fit(module,new,args.out/'corrected')
        write_json(common.compare(old_fit,new_fit),args.out/'comparison.json')
        write_json(dict(old=str(old_fit),corrected=str(new_fit),
                        biological_specimens=1,independent_validation=False),args.out/'report.json')
        common.write_inventory(args.out)
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
