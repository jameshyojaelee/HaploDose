"""Exhaustive independent cell-pair enumeration for DNA stability summaries."""
import argparse
from collections import Counter
import importlib.util
import itertools
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    helper=Path(__file__).with_name('andor_dna_contrast_stability.py')
    spec=importlib.util.spec_from_file_location('contrast',helper)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    run=start_run(args.out,dict(script=Path(__file__),helper=helper,design=args.task/'CONTRAST_DESIGN.md'),{},
                  dict(scope='Exhaustive integer fixtures, biological n0'),mode='andor_dna_contrast_preflight')
    try:
        anchor=['g0','g0','g1','g1','g1'];tested=0
        for candidate in itertools.product(('c0','c1','excluded','unassigned'),repeat=5):
            expected=Counter()
            for i,j in itertools.combinations(range(5),2):
                left,right=sorted((anchor[i],anchor[j]));key=(left,right)
                expected[key+('fixed_reference_pairs',)]+=1
                valid=candidate[i] not in ('excluded','unassigned') and candidate[j] not in ('excluded','unassigned')
                expected[key+('both_assigned_pairs',)]+=int(valid)
                expected[key+('coassigned_pairs',)]+=int(valid and candidate[i]==candidate[j])
                expected[key+('separated_pairs',)]+=int(valid and candidate[i]!=candidate[j])
                expected[key+('unresolved_pairs',)]+=int(not valid)
            actual=module.pair_counts(anchor,candidate)
            for row in actual:
                for field in ('fixed_reference_pairs','both_assigned_pairs','coassigned_pairs','separated_pairs','unresolved_pairs'):
                    assert row[field]==expected[(row['anchor_left'],row['anchor_right'],field)]
            renamed=tuple({'c0':'z9','c1':'z1'}.get(value,value) for value in candidate)
            assert actual==module.pair_counts(anchor,renamed)
            assert actual==module.pair_counts(anchor[::-1],candidate[::-1]);tested+=1
        source=args.out/'invalid_cells.tsv'
        reference_unresolved=module.pair_counts(['excluded','g0','g0','unassigned'],['c0']*4)
        assert len(reference_unresolved)==1 and reference_unresolved[0]['fixed_reference_pairs']==1
        assert reference_unresolved[0]['coassigned_pairs']==1
        source.write_text('rg\tclone\n001\tc0\nNA\tc1\n')
        assert module.read_groups(source).index.tolist()==['001','NA']
        for content in ('rg\tclone\n001\tc0\n001\tc1\n','rg\tclone\n\tc0\n','rg\tclone\n001\t\n'):
            source.write_text(content)
            try:
                module.read_groups(source)
            except AssertionError:
                pass
            else:
                raise AssertionError('Ambiguous cell identity/group accepted')
        write_json(dict(exhaustive_candidate_assignments=tested,independent_pair_enumeration_exact=True,
                        label_and_row_invariance_exact=True,invalid_identity_or_group_controls=3,
                        biological_specimens=0),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
