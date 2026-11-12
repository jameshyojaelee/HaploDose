"""Postfit membership/exposure factorial attribution; no oracle inputs to inference."""
import argparse, importlib.util, json, sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.optimize import minimize, brentq
from haplodose.simulate import simulate
from haplodose.provenance import start_run, finish_run


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--trace-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    path=Path(__file__).with_name('membership_latent_fit.py')
    spec=importlib.util.spec_from_file_location('conditional_reference',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    rec=start_run(a.out,{'runner':__file__,'likelihood':path},{},
        {'scope':'Postfit oracle attribution; not deployable performance or calibrated confidence'},mode='attribution')
    results=[]
    try:
        for seed in (20261108,20261109,20261110):
            capture=[]
            def observe(frame,event,result):
                if event=='return' and frame.f_code.co_name=='draw' and frame.f_code.co_filename==simulate.__code__.co_filename:
                    loc=frame.f_locals
                    if loc.get('cols',{}).get('specimen')=='S1':
                        capture.append((loc['cols']['truth_clone'],loc['L'][:,None]*loc['share']))
            sys.setprofile(observe)
            try:X,genes,meta,*_=simulate(seed=seed)
            finally:sys.setprofile(None)
            source=a.panel_root/f'seed_{seed}'/'data/counts.h5'
            with h5py.File(source) as h:
                for name,values in [('data',X.data),('indices',X.indices),('indptr',X.indptr)]:
                    if not np.array_equal(values,h['X/'+name][:]):raise ValueError('Frozen counts do not reproduce')
                if not np.array_equal(meta.cell_id.to_numpy(),h['cell_ids'].asstr()[:]):raise ValueError('Cell identity mismatch')
            if [x[0] for x in capture]!=['normal','clone1','clone2']:raise ValueError('Unexpected population order')
            mu=np.vstack([x[1] for x in capture]);s1=np.flatnonzero(meta.specimen.eq('S1'))
            root=a.trace_root/f'trace_{seed}'
            mapping=pd.read_csv(root/'fit/genes.tsv.gz',sep='\t');mapped=mapping[mapping.bin>=0]
            bmean=np.zeros((len(mu),int(mapped.bin.max())+1))
            for col,binid in zip(mapped.column.to_numpy(int),mapped.bin.to_numpy(int)):
                bmean[:,binid]+=mu[:,col]
            for event in json.loads((root/'trace.json').read_text()):
                if event['stage']!='event_fit':continue
                z=np.load(root/event['arrays']);global_rows=z['rows'][z['cells']]
                local=np.searchsorted(s1,global_rows)
                if not np.array_equal(s1[local],global_rows):raise ValueError('Specimen alignment failure')
                y=(meta.iloc[global_rows].truth_clone.to_numpy()=='clone2').astype(float)
                delta=.6*np.log(1.5)
                # Every saved event bin is in the implanted first half of chr12.
                selected=mapped[mapped.bin.isin(z['cols'])]
                positions=pd.read_csv(a.panel_root/f'seed_{seed}'/'data/positions.tsv',sep='\t')
                if not np.array_equal(positions.gene.to_numpy(),np.asarray(genes)):
                    raise ValueError('Position/gene order mismatch')
                positions=positions.iloc[selected.column.to_numpy(int)]
                if not (positions.chromosome.eq('chr12').all() and (positions.start<=10000000).all()):
                    raise ValueError('Oracle effect is not constant over selected event bins')
                Eoracle=.5*bmean[local][:,z['cols']]/np.exp(delta*y[:,None])
                for ename,E in [('fitted',z['E']),('oracle',Eoracle)]:
                    for qname,q in [('soft',z['q']),('hard',(z['q']>.5).astype(float)),('truth',y)]:
                        O,r=z['Ytest'],z['r'];fun=lambda x:m.objective(x,O,E,r,q)[:2]
                        initial=z['parameters'][:2]
                        opt=minimize(fun,initial,jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':1000})
                        ref=minimize(lambda x:fun(x)[0],initial,method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-10,'maxiter':2000})
                        value,grad,H,_=m.objective(opt.x,O,E,r,q);eig=np.linalg.eigvalsh(H);direction=np.array([-1.,1.])
                        roots=None
                        if qname!='soft':
                            roots=[]
                            for group in (0,1):
                                take=q==group
                                rg=np.broadcast_to(r,O.shape)[take]
                                def score(x):
                                    mean=E[take]*np.exp(x)
                                    return float((rg*(O[take]-mean)/(rg+mean)).sum())
                                lo,hi=-1.,1.
                                while score(lo)<=0:lo*=2
                                while score(hi)>=0:hi*=2
                                roots.append(brentq(score,lo,hi,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps))
                        row=dict(seed=seed,exposure=ename,membership=qname,known_log_effect=delta,
                            fitted_log_effect=float(direction@opt.x),conditional_SE=float(np.sqrt(direction@np.linalg.solve(H,direction))) if min(eig)>0 else None,
                            parameters=opt.x.tolist(),gradient=grad.tolist(),curvature_eigenvalues=eig.tolist(),
                            optimizer_success=bool(opt.success),optimizer_message=str(opt.message),
                            independent_optimizer_success=bool(ref.success),independent_objective_discrepancy=float(value-ref.fun),
                            independent_scalar_roots=roots,oracle_is_deployable=False,
                            scope='Conditional working-NB likelihood; fitted q, nuisance, aggregation, thinning misspecification and selection uncertainty excluded')
                        results.append(row);print(json.dumps(row),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise

if __name__=='__main__':main()
