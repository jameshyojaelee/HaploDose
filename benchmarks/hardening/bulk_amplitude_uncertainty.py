"""Indexed joint reference/calibration/patient bootstrap for frozen bulk amplitudes."""
import argparse
from collections import Counter
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

SEED=20261109
TOTAL=1000000
CANDIDATES=('hd_path_mean','hd_path_precision_mean')
ENDPOINTS=('pearson','recall_loss','fpr_loss','precision_loss')


def load(name):
    spec=importlib.util.spec_from_file_location('amplitude_uncertainty_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


class Calibration:
    """Exact equal-patient quantile using precomputed integer cumulative masses."""
    def __init__(self,scores,neutral,columns):
        scores=np.asarray(scores,float);neutral=np.asarray(neutral,bool)
        if scores.ndim!=2 or scores.shape!=neutral.shape:raise ValueError('Unaligned calibration arrays')
        self.columns=np.asarray(columns,int)
        if len(np.unique(self.columns))!=len(self.columns):raise ValueError('Unique calibration population required')
        values=[np.sort(scores[:,j][neutral[:,j]&np.isfinite(scores[:,j])]) for j in self.columns]
        self.present=np.array([len(v)>0 for v in values])
        self.denominator=math.lcm(*(len(v) for v in values if len(v)))
        self.support=np.unique(np.concatenate(values)) if values else np.empty(0)
        self.mass=np.empty((len(self.support),len(values)),object)
        for j,v in enumerate(values):
            weight=self.denominator//len(v) if len(v) else 0
            self.mass[:,j]=[weight*int(n) for n in np.searchsorted(v,self.support,side='right')]
        self.position={int(j):k for k,j in enumerate(self.columns)}

    def threshold(self,chosen):
        counts=Counter(int(j) for j in chosen)
        if any(j not in self.position for j in counts):raise ValueError('Calibration identity outside frozen population')
        entries=[(self.position[j],n) for j,n in counts.items() if self.present[self.position[j]]]
        if not entries:return np.nan
        # Fixed decimal 5% = 1/20. Absent patients contribute no calibration mass.
        target=self.denominator*sum(n for _,n in entries)
        lo,hi=0,len(self.support)-1
        while lo<hi:
            mid=(lo+hi)//2
            mass=sum(self.mass[mid,j]*n for j,n in entries)
            if 20*mass>=target:hi=mid
            else:lo=mid+1
        return float(self.support[lo])


class LossCurve:
    """Strict-threshold recall/FPR/precision, preserving distinct denominators."""
    def __init__(self,scores,loss,neutral,gain,common):
        scores=np.asarray(scores,float)
        loss,neutral,gain,common=(np.asarray(x,bool) for x in (loss,neutral,gain,common))
        if any(x.shape!=scores.shape for x in (loss,neutral,gain,common)):raise ValueError('Unaligned scoring arrays')
        if np.any(loss&(neutral|gain)):raise ValueError('Overlapping positive and negative truth')
        self.denominators=np.stack([(loss&common).sum(0),(neutral&common).sum(0)])
        self.support=np.unique(scores[common&~np.isnan(scores)])
        finite=np.isfinite(scores)
        masks=(loss&common,neutral&common,loss&common&finite,(loss|neutral|gain)&common&finite)
        self.counts=np.empty((len(self.support)+1,4,scores.shape[1]),np.int64)
        for k,mask in enumerate(masks):
            for j in range(scores.shape[1]):
                v=np.sort(scores[:,j][mask[:,j]&~np.isnan(scores[:,j])])
                self.counts[:-1,k,j]=np.searchsorted(v,self.support,side='left')
                self.counts[-1,k,j]=len(v)

    def at(self,threshold):
        out=np.full((self.counts.shape[2],3),np.nan)
        if np.isfinite(threshold):
            counts=self.counts[np.searchsorted(self.support,threshold,side='left')]
            for k in range(2):np.divide(counts[k],self.denominators[k],out=out[:,k],where=self.denominators[k]>0)
            np.divide(counts[2],counts[3],out=out[:,2],where=counts[3]>0)
        return out


def definitions(methods):
    controls=[m for m in methods if m not in CANDIDATES]
    assert len(methods)==12 and len(controls)==10
    return [dict(method=m,control=c,endpoint=e,direction=-1 if e=='fpr_loss' else 1)
            for m in CANDIDATES for c in controls for e in ENDPOINTS]


def effects(scores,chosen,contrasts):
    values=np.array([(scores[r['method']][:,ENDPOINTS.index(r['endpoint'])]-
                      scores[r['control']][:,ENDPOINTS.index(r['endpoint'])])[chosen]*r['direction'] for r in contrasts])
    finite=np.isfinite(values);n=finite.sum(1);result=np.full(len(contrasts),np.nan)
    np.divide(np.where(finite,values,0).sum(1),n,out=result,where=n>0)
    return result,n


def read_bank(root):
    bank=root/'dev/bulk_amplitude_reference_20261109';windows=[];correlations=[];inputs={}
    report=json.loads((bank/'summary_v1/report.json').read_text());assert report['all_nine_independent_audits_passed']
    inputs['bank_audit']=bank/'summary_v1/report.json'
    methods=json.loads((bank/'refit_000_v1/report.json').read_text())['methods']
    for i in range(9):
        directory=bank/f'refit_{i:03d}_v1';inputs[f'windows_{i}']=directory/'window_scores.npz';inputs[f'scores_{i}']=directory/'patient_scores.tsv'
        inputs[f'audit_{i}']=bank/f'audit_{i:03d}_v1/report.json'
        assert json.loads(inputs[f'audit_{i}'].read_text())['all_checks_passed']
        with np.load(inputs[f'windows_{i}']) as z:w={k:z[k] for k in z.files}
        if i:
            for key in ('samples','common','loss','neutral','gain','dna','split'):
                assert np.array_equal(w[key],windows[0][key],equal_nan=True) if w[key].dtype.kind not in 'US' else np.array_equal(w[key],windows[0][key])
        table=pd.read_csv(inputs[f'scores_{i}'],sep='\t',float_precision='round_trip').set_index(['sample','method'])
        assert table.index.is_unique
        correlations.append({m:table.xs(m,level='method').reindex(w['samples']).pearson.to_numpy(float) for m in methods})
        windows.append(w)
    return windows,correlations,methods,inputs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--start',type=int,required=True);p.add_argument('--stop',type=int,required=True)
    p.add_argument('--preflight-report',type=Path)
    a=p.parse_args();assert 0<=a.start<a.stop<=TOTAL;a.out.mkdir(exist_ok=False)
    root=a.task.parents[1];bank,corr,methods,inputs=read_bank(root)
    if a.stop>1000:
        assert a.preflight_report and json.loads(a.preflight_report.read_text())['all_checks_passed']
        inputs['preflight']=a.preflight_report
    inputs.update(script=Path(__file__),design=a.task/'DESIGN.md',profile=Path(__file__).with_name('profile.py'))
    rec=start_run(a.out,inputs,dict(start=a.start,stop=a.stop,total=TOTAL,seed=SEED,scientific_family=80),
        dict(scope='Approximate joint bootstrap conditional on eight frozen donor-reference draws; not validated coverage or independent evaluation'),mode='bulk_amplitude_uncertainty')
    try:
        ci=np.flatnonzero(bank[0]['split']=='calibration');ei=np.flatnonzero(bank[0]['split']=='evaluation');assert len(ci)==len(ei)==24
        contrasts=definitions(methods);assert len(contrasts)==80
        calibrations=[];curves=[]
        for w in bank:
            calibrations.append({m:Calibration(w[m],w['neutral']&w['common'],ci) for m in methods})
            curves.append({m:LossCurve(w[m],w['loss'],w['neutral'],w['gain'],w['common']) for m in methods})
        def scored(ri,cal):
            thresholds={m:calibrations[ri][m].threshold(cal) for m in methods}
            return {m:np.column_stack([corr[ri][m],curves[ri][m].at(thresholds[m])]) for m in methods},thresholds
        point,_=scored(0,ci);point_effects,point_n=effects(point,ei,contrasts)
        (a.out/'definition.json').write_text(json.dumps(dict(methods=methods,contrasts=contrasts,
            point_effects=point_effects.tolist(),point_paired_patients=point_n.tolist(),
            effect_orientation='Positive favours candidate; FPR difference sign reversed',seed=SEED,total=TOTAL),indent=2)+'\n')
        size=a.stop-a.start
        draws=np.lib.format.open_memmap(a.out/'draws.npy',mode='w+',dtype=float,shape=(size,80))
        paired=np.lib.format.open_memmap(a.out/'paired.npy',mode='w+',dtype=np.uint8,shape=(size,80))
        indices=np.lib.format.open_memmap(a.out/'indices.npy',mode='w+',dtype=np.int16,shape=(size,49))
        thresholds=np.lib.format.open_memmap(a.out/'thresholds.npy',mode='w+',dtype=float,shape=(size,12))
        rng=np.random.default_rng(SEED)
        for k in range(a.stop):
            ri=int(rng.integers(1,9));cal=rng.choice(ci,24,replace=True);chosen=rng.choice(ei,24,replace=True)
            if k<a.start:continue
            scores,t=scored(ri,cal);row=k-a.start
            draws[row],paired[row]=effects(scores,chosen,contrasts)
            indices[row]=np.r_[ri,cal,chosen];thresholds[row]=[t[m] for m in methods]
            if (k+1)%1000==0:print(json.dumps(dict(completed=k+1,start=a.start,stop=a.stop)),flush=True)
        for x in (draws,paired,indices,thresholds):x.flush()
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        (a.out/'report.json').write_text(json.dumps(dict(start=a.start,stop=a.stop,scientific_family=80,
            method_promoted=False,coverage_validated=False),indent=2)+'\n')
        finish_run(a.out,rec,backend='CPU exact integer calibration and strict-threshold score lookup')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
