"""Two-factor RNA effect diagnostic; median outputs are not count observations."""
import math
import numpy as np
from scipy import sparse
from haplodose.profile import log_ratio

NAMES=('pooled_sum','equal_sum','pooled_median','equal_median')


def effects(counts,bin_ids,normal):
    """Genes x samples; equal-donor reference is arithmetic mean of normal rates.

    Sum bins use the existing half-count offset per bin. Median bins use the
    existing bulk median-mode half-count offset per gene. This change of the
    elementary observation and its offset is part of the aggregation intervention.
    """
    G=np.asarray(counts,float);b=np.asarray(bin_ids,int);normal=np.asarray(normal,bool)
    if G.ndim!=2 or b.shape!=(len(G),) or normal.shape!=(G.shape[1],):raise ValueError('Misaligned inputs')
    if not np.isfinite(G).all() or np.any(G<0) or not np.equal(G,np.floor(G)).all():raise ValueError('Raw counts required')
    if np.any(b<0) or not np.array_equal(np.unique(b),np.arange(b.max()+1)):raise ValueError('Nonconsecutive bin IDs')
    if normal.sum()<2:raise ValueError('At least two reference donors required')
    A=sparse.csr_matrix((np.ones(len(G)),(np.arange(len(G)),b)),shape=(len(G),int(b.max()+1)))
    Y=np.asarray((A.T@G).T);L=Y.sum(axis=1)
    np.testing.assert_array_equal(L,G.sum(axis=0))
    if np.any(L<=0) or np.any(G[:,normal].sum(axis=1)<=0):raise ValueError('Undefined reference rates')
    gene_rates=dict(pooled=G[:,normal].sum(axis=1)/L[normal].sum(),
                    equal=(G[:,normal]/L[normal]).mean(axis=1))
    # Use direct bin arithmetic for the existing control to preserve its exact
    # evaluation order; aggregation of gene rates is mathematically equivalent.
    bin_rates=dict(pooled=Y[normal].sum(axis=0)/Y[normal].sum(),
                   equal=(Y[normal]/L[normal,None]).mean(axis=0))
    output={}
    for ref,mu in gene_rates.items():
        aggregate=np.asarray(A.T@mu).ravel()
        allowance=2*len(G)*np.finfo(float).eps*np.maximum(aggregate,bin_rates[ref])
        if np.any(np.abs(aggregate-bin_rates[ref])>allowance):raise ArithmeticError('Reference aggregation identity failed')
        summed=log_ratio(Y,L[:,None]*bin_rates[ref][None])
        gene_lr=np.log((G+.5)/(mu[:,None]*L[None]+.5))
        median=np.column_stack([np.median(gene_lr[b==j],axis=0) for j in range(Y.shape[1])])
        for kind,values in [('sum',summed),('median',median)]:
            output[ref+'_'+kind]=values-np.median(values,axis=1,keepdims=True)
    return output,dict(depth=L,reference_gene_rates=gene_rates,reference_bin_rates=bin_rates)


def verify():
    G=np.array([[90,1,40],[10,9,20],[20,10,30],[10,20,10]],float)
    b=np.array([0,0,1,1]);normal=np.array([True,True,False]);actual,meta=effects(G,b,normal)
    # Independent scalar calculation of both reference constructions and effects.
    L=[math.fsum(G[:,s]) for s in range(3)]
    scalar={}
    for ref in ('pooled','equal'):
        mu=[math.fsum(G[g,:2])/math.fsum(L[:2]) if ref=='pooled' else
            math.fsum(G[g,s]/L[s] for s in range(2))/2 for g in range(4)]
        for kind in ('sum','median'):
            matrix=[]
            for s in range(3):
                row=[]
                for j in range(2):
                    genes=np.flatnonzero(b==j)
                    value=(math.log((math.fsum(G[genes,s])+.5)/(L[s]*math.fsum(mu[g] for g in genes)+.5)) if kind=='sum' else
                        float(np.median([math.log((G[g,s]+.5)/(L[s]*mu[g]+.5)) for g in genes])))
                    row.append(value)
                matrix.append(np.asarray(row)-np.median(row))
            scalar[ref+'_'+kind]=np.asarray(matrix)
    allowance=32*np.finfo(float).eps*max(1.,max(float(np.max(np.abs(v))) for v in scalar.values()))
    for key in NAMES:np.testing.assert_allclose(actual[key],scalar[key],rtol=0,atol=allowance)
    perm=np.array([2,0,3,1]);reordered,_=effects(G[perm],b[perm],normal)
    for key in NAMES:np.testing.assert_allclose(actual[key],reordered[key],rtol=0,atol=allowance)
    doubled=G.copy();doubled[:,0]*=2;_,other=effects(doubled,b,normal)
    np.testing.assert_allclose(meta['reference_gene_rates']['equal'],other['reference_gene_rates']['equal'],rtol=0,atol=allowance)
    if np.allclose(meta['reference_gene_rates']['pooled'],other['reference_gene_rates']['pooled']):
        raise AssertionError('Fixture failed to distinguish depth weighting')
    return dict(independent_scalar_reference=True,gene_permutation=True,equal_reference_depth_invariance=True,
                roundoff_allowance=allowance)
