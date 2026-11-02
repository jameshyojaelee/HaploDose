"""Cohort table for the U-CAN Visium HD vs FACETS benchmark. Usage: summarise.py ROOT
Patient is the unit (one section per patient). Adds IGF2 TD carrier status from the CRC WGS analysis table."""
import json
import sys
from pathlib import Path

import pandas as pd

TD = Path('/gpfs/commons/home/jameslee/CRC/results/discovery_a19_sv_outcome_base/analysis_table.tsv')

root = Path(sys.argv[1])
sub = sys.argv[2] if len(sys.argv) > 2 else 'dna'
rows = [json.loads(p.read_text()) for p in sorted(root.glob(f'*/{sub}/summary.json'))]
t = pd.DataFrame(rows).drop(columns=['vcf', 'thresholds'])
td = pd.read_csv(TD, sep='\t', usecols=['patient_id', 'td_carrier'])
td['sample'] = td.patient_id.str.extract(r'(U\d{4})', expand=False)
t = t.merge(td[['sample', 'td_carrier']].drop_duplicates('sample'), on='sample', how='left')
t['wgd_like'] = t.ploidy >= 2.7          # CRC R7-06 definition
t['altered_recall'] = t.altered_arms_same_direction / t.dna_altered_arms
t['neutral_false_rate'] = t.neutral_arms_called_altered / t.neutral_arms
t.to_csv(root / f'cohort_summary_{sub}.tsv', sep='\t', index=False)
cols = ['sample', 'td_carrier', 'ploidy', 'purity', 'tumour_clones', 'rho_log_ratio', 'dna_altered_arms',
        'altered_arms_same_direction', 'neutral_arms_called_altered', 'subclonal_segments', 'subclonal_segments_split']
print(t[cols].to_string(index=False, float_format=lambda x: f'{x:.2f}'))
tot = t[['dna_altered_arms', 'altered_arms_same_direction', 'neutral_arms', 'neutral_arms_called_altered',
         'subclonal_segments', 'subclonal_segments_split']].sum()
print(f"\nn patients {len(t)} | median rho {t.rho_log_ratio.median():.2f} | altered arms right direction "
      f"{tot.altered_arms_same_direction}/{tot.dna_altered_arms} | neutral arms called altered "
      f"{tot.neutral_arms_called_altered}/{tot.neutral_arms} | subclonal DNA segments split "
      f"{tot.subclonal_segments_split}/{tot.subclonal_segments}")
t['subclonal_chance_hits'] = t.subclonal_split_fraction_expected_by_chance * t.subclonal_segments
print(f"subclonal splits observed {int(t.subclonal_segments_split.sum())} vs expected by chance "
      f"{t.subclonal_chance_hits.sum():.1f} of {int(t.subclonal_segments.sum())}")
print('median rho by WGD-like:', t.groupby('wgd_like').rho_log_ratio.median().round(2).to_dict())
