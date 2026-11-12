"""Paired old/exact threshold replay on all frozen factorial profiles."""
import argparse,importlib.util,json
from fractions import Fraction
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];factor=root/'dev/gene_factorial_20261108';panel=root/'dev/hardening_20261106_2119/profile_v6'
 new=load('exact_profile',Path(__file__).with_name('profile.py'));oldpath=a.task/'source_v1/benchmarks/hardening/profile.py';old=load('old_profile',oldpath)
 manifest=json.loads((factor/'panel_v1/run_manifest.json').read_text());assert manifest['status']=='complete'
 inputs={k:Path(v['path']) for k,v in manifest['input']['files'].items()};inputs.update(script=Path(__file__),new_scorer=Path(__file__).with_name('profile.py'),old_scorer=oldpath)
 for i in range(9):inputs[f'factor{i}']=factor/f'panel_v1/refit_{i:03d}.npz'
 rec=start_run(a.out,inputs,dict(fpr='1/20',calls='strict score<threshold'),dict(scope='Corrected development loss scoring; no model change'),mode='loss_quantile_replay')
 try:
  genes=(panel/'gene_universe.txt').read_text().splitlines();tab,seg=new.subarm.absolute();truth=(new.subarm.hg19_genes(),tab,seg,*new.subarm.hg19_arms());d=new.profiles(panel/'raw_refit_000',panel/'raw_refit_000/haplodose',truth,genes)
  with np.load(panel/'raw_support_and_truth.npz',allow_pickle=True) as z:common=z['common'];samples=z['samples']
  np.testing.assert_array_equal(samples,d['tum']);meta=pd.read_csv(panel/'split_manifest.tsv',sep='\t').set_index('sample').loc[samples];ci=np.flatnonzero(meta.split.eq('calibration'));assert len(ci)==24
  comparisons=[];rows=[]
  for i in range(9):
   with np.load(inputs[f'factor{i}'],allow_pickle=True) as z:
    np.testing.assert_array_equal(z['samples'],samples)
    for method in [n for n in z.files if n not in ('samples','arms')]:
     scores=z[method];neutral=d['NEU']&common
     before=old.loss_threshold(scores,neutral,ci);after=new.loss_threshold(scores,neutral,ci)
     below=Fraction(0);at=Fraction(0);patients=0
     for j in ci:
      v=scores[:,j][neutral[:,j]&np.isfinite(scores[:,j])]
      if len(v):
       patients+=1;below+=Fraction(int(np.sum(v<after)),len(v));at+=Fraction(int(np.sum(v<=after)),len(v))
     assert patients and below<Fraction(patients,20)<=at
     one=old.score_loss(scores,d['FL'],d['NEU'],common,before,negative=d['FG']|d['NEU']);two=new.score_loss(scores,d['FL'],d['NEU'],common,after,negative=d['FG']|d['NEU'])
     if before==after:pd.testing.assert_frame_equal(pd.DataFrame(one),pd.DataFrame(two),check_exact=True)
     comparisons.append(dict(refit=i,method=method,before=before,after=after,changed=before!=after,calibration_patients=patients,exact_mass_below=str(below/patients),exact_mass_at=str(at/patients)))
     for mode,items in [('historical',one),('exact',two)]:
      for sample,row in zip(samples,items):rows.append(dict(refit=i,method=method,scorer=mode,sample=sample,split=meta.loc[sample,'split'],**row))
  pd.DataFrame(comparisons).to_csv(a.out/'thresholds.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
  for item in rec['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  (a.out/'report.json').write_text(json.dumps(dict(thresholds_checked=len(comparisons),thresholds_changed=sum(v['changed'] for v in comparisons),patient_score_rows=len(rows),all_exact_CDF_checks_pass=True,scientific_target_changed=False,scope='Point replay only; bootstrap correction remains required'),indent=2))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
