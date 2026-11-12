"""Prespecified fresh full-workflow comparison of fixed research gene-veto arms."""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.simulate import simulate


CONDITIONS = ('matched', 'state', 'two_clone')
MASTER = 20261109


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def generate(rep, phase, condition):
    source_seed = MASTER*10000 + phase*1000 + rep
    capture = []
    def observe(frame, event, value):
        if event == 'return' and frame.f_code.co_name == 'draw' and frame.f_code.co_filename == simulate.__code__.co_filename:
            loc = frame.f_locals
            if loc['cols']['specimen'] == 'S1':
                capture.append((loc['L'].copy(), loc['share'].copy()))
    assert sys.getprofile() is None
    sys.setprofile(observe)
    try:
        X, genes, metadata, positions, _, _ = simulate(seed=source_seed,
            specimens=[dict(name='S1',n_cells=600,normal_frac=0.,clones=[dict(frac=1.,events=[
                ('chr3',0.,1.,3),('chr8',.5,1.,1)])])])
    finally:
        sys.setprofile(None)
    assert len(capture) == 1
    exposure, share = capture[0]; base = exposure[:,None]*share
    query = np.flatnonzero(metadata.specimen.eq('S1'))
    assert len(query) == 600 and base.shape == (600,len(genes))
    delta = .6*np.log(1.5); sigma = delta*np.sqrt(.6*.4)
    rng_state = np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,1]))
    state = rng_state.normal(0.,sigma,600)
    rng_labels = np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,2]))
    labels = np.zeros(600,int); labels[rng_labels.permutation(600)[:240]] = 1
    region = positions.chromosome.eq('chr12').to_numpy() & (positions.start.to_numpy() <= 10000000)
    effect = np.zeros(600) if condition == 'matched' else state.copy()
    if condition == 'two_clone':
        effect += labels*delta
    mean = base*np.exp(effect[:,None]*region[None])
    mean *= (exposure/mean.sum(axis=1))[:,None]
    gamma = np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,3])).gamma(10.,.1,mean.shape)
    counts = np.random.default_rng(np.random.SeedSequence([MASTER,phase,rep,4])).poisson(mean*gamma).astype(np.int32)
    # Reference counts remain the original source measurements in every condition.
    keep = np.flatnonzero(~metadata.specimen.eq('S1'))
    assert np.array_equal(np.r_[keep,query],np.arange(len(metadata)))
    X = sparse.vstack([X[keep],sparse.csr_matrix(counts)]).tocsr()
    truth = labels if condition == 'two_clone' else np.zeros(600,int)
    metadata.loc[query,'truth_clone'] = np.where(truth == 0,'clone1','clone2')
    arrays = dict(base_mean=base,mean=mean,exposure=exposure,state=state,latent_labels=labels,
        truth=truth,region=region,query_rows=query,delta=np.array(delta),sigma=np.array(sigma))
    return source_seed, X, genes, metadata, positions, arrays


def write_input(directory, X, genes, metadata, positions):
    directory.mkdir()
    with h5py.File(directory/'counts.h5','w') as h:
        for name,array in (('data',X.data),('indices',X.indices),('indptr',X.indptr)):
            h.create_dataset('X/'+name,data=array,compression='gzip')
        h.create_dataset('genes',data=np.array(genes,dtype=object),dtype=h5py.string_dtype())
        h.create_dataset('cell_ids',data=metadata.cell_id.to_numpy(dtype=object),dtype=h5py.string_dtype())
    metadata.to_csv(directory/'truth_cells.tsv',sep='\t',index=False)
    metadata.drop(columns=[c for c in metadata if c.startswith('truth_')]).to_csv(directory/'model_cells.tsv',sep='\t',index=False)
    positions.to_csv(directory/'positions.tsv',sep='\t',index=False)


def score(folder, metadata, condition):
    cells = pd.read_csv(folder/'cells.tsv.gz',sep='\t')
    truth = metadata.loc[metadata.specimen.eq('S1'),['cell_id','truth_clone']]
    cells = cells[cells.specimen.eq('S1')].merge(truth,on='cell_id',validate='one_to_one')
    assert len(cells) == 600 and not cells.truth_clone.isna().any()
    profiles = pd.read_csv(folder/'clones.tsv',sep='\t')
    profiles = profiles[profiles.specimen.eq('S1') & ~profiles.normal_like]
    names = sorted(profiles.clone)
    differences = pd.read_csv(folder/'differences.tsv',sep='\t')
    differences = differences[differences.specimen.eq('S1') & differences.clone_a.isin(names) & differences.clone_b.isin(names)]
    chromosomes = sorted(set(differences.chromosome))
    populations = sorted(set(cells.truth_clone))
    table = np.array([[sum((cells.truth_clone == p)&(cells.clone == g)&cells.malignant_call)
        for g in names] for p in populations],dtype=int)
    correct = int(table[linear_sum_assignment(-table)].sum()) if table.size else 0
    return dict(condition=condition,reported_groups=len(names),null_error=bool(len(names)>1),
        exact_chr12_recovery=bool(len(names)==2 and chromosomes == ['chr12']),
        any_two_groups=bool(len(names)>=2),difference_chromosomes=chromosomes,
        assignment_correct=correct,assignment_n=600,tumour_calls=int(cells.malignant_call.sum()),
        assignment_table=table.tolist(),populations=populations,groups=names)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--phase',type=int,choices=(0,1),required=True)
    p.add_argument('--condition',choices=CONDITIONS+('all',),required=True)
    p.add_argument('--start',type=int,required=True);p.add_argument('--stop',type=int,required=True)
    a = p.parse_args();a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    hook = load('fresh_gene_hook',root/'gene_thinning_hook.py')
    stage = load('fresh_stage',root/'membership_weighted_stages.py')
    veto = load('fresh_veto',root/'gene_veto_weights.py')
    inputs = dict(script=Path(__file__),hook=root/'gene_thinning_hook.py',stage=root/'membership_weighted_stages.py',
        veto=root/'gene_veto_weights.py',design=a.task/'DESIGN.md',
        candidate_gate=a.task.parent/'gene_veto_aligned_20261109/audit_v2/report.json')
    params = clones.CloneParams(split_mode='thin')
    rec = start_run(a.out,inputs,dict(phase=a.phase,start=a.start,stop=a.stop,condition=a.condition,params=asdict(params)),
        dict(scope='Fresh conditional known-noise gene experiment, biological n=0'),mode='gene_veto_fresh')
    results = []
    try:
        assert json.loads(inputs['candidate_gate'].read_text())['all_checks_passed']
        conditions = CONDITIONS if a.condition == 'all' else (a.condition,)
        for rep in range(a.start,a.stop):
            for condition in conditions:
                seed,X,genes,metadata,positions,arrays = generate(rep,a.phase,condition)
                case = a.out/condition/f'seed_{seed}';case.mkdir(parents=True)
                write_input(case/'input',X,genes,metadata,positions)
                np.savez_compressed(case/'generation.npz',**arrays)
                for arm in ('comparator','candidate'):
                    out = case/arm;out.mkdir();folds=out/'folds';folds.mkdir();vdir=out/'veto';vdir.mkdir()
                    gene_trace,stage_trace,veto_trace = [],[],[]
                    with hook.gene_native(gene_trace,folds,case_seeds=(seed,)), stage.weighted_stages('both',stage_trace), \
                            veto.observe_weights(veto_trace,vdir,apply=arm == 'candidate'):
                        pipeline.run(case/'input/counts.h5',case/'input/model_cells.tsv',case/'input/positions.tsv',out/'fit',
                            params=params,chunk=97,counts_kind='raw',genome_build='hg38',coordinate_system='one-based-closed',gene_namespace='symbol')
                    write_json(dict(gene=gene_trace,stage=stage_trace,veto=veto_trace),out/'trace.json')
                    row=score(out/'fit',metadata,condition)
                    row.update(phase=a.phase,rep=rep,seed=seed,arm=arm,case=str(case),
                        fit_manifest=file_record(out/'fit/run_manifest.json'))
                    results.append(row);write_json(results,a.out/'results.json')
                    print(json.dumps(row),flush=True)
        finish_run(a.out,rec,backend='CPU whole gene-count workflow; fixed comparator/candidate')
    except BaseException as error:
        finish_run(a.out,rec,error=error);raise


if __name__ == '__main__':
    main()
