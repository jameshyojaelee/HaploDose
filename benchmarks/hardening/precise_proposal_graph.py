"""Research numerical correction: float64 PCA and identity-based self exclusion.

State grids, neighbourhood size, randomized PCA settings and Leiden parameters
are unchanged. No biological reliability claim follows from graph reproducibility.
"""
from contextlib import contextmanager
import random
import igraph as ig
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from haplodose import clones


def nonself_neighbors(X,k):
    candidates=NearestNeighbors(n_neighbors=k+1).fit(X).kneighbors(X,return_distance=False)
    # With duplicate features, a query's own row need not be the first result
    # or even appear in this bounded candidate list. Remove by identity only.
    return np.array([row[row!=i][:k] for i,row in enumerate(candidates)])


def precise_leiden_groups(F,params):
    n=F.shape[0]
    if n<2*params.min_cells:return np.zeros(n,dtype=int)
    X=PCA(n_components=int(min(params.n_pcs,n-1,F.shape[1]-1)),svd_solver='randomized',
          random_state=params.seed).fit_transform(np.asarray(F,np.float64))
    k=int(min(params.n_neighbors,n-1));idx=nonself_neighbors(X,k)
    edges=np.column_stack([np.repeat(np.arange(n),k),idx.ravel()])
    g=ig.Graph(n=n,edges=edges,directed=False).simplify()
    ig.set_random_number_generator(random.Random(params.seed))
    part=g.community_leiden(objective_function='modularity',resolution=params.leiden_resolution,n_iterations=2)
    return clones.consolidate(np.asarray(part.membership),X,params)


@contextmanager
def precise_graph():
    old=clones.leiden_groups;clones.leiden_groups=precise_leiden_groups
    try:yield
    finally:clones.leiden_groups=old


def checks():
    duplicate=np.zeros((6,2));idx=nonself_neighbors(duplicate,2)
    if idx.shape!=(6,2) or any(i in row or len(set(row))!=2 for i,row in enumerate(idx)):
        raise AssertionError('Duplicate-feature self exclusion failed')
    X=np.random.default_rng(20261108).normal(size=(31,7));idx=nonself_neighbors(X,5)
    distances=cdist(X,X,metric='sqeuclidean');np.fill_diagonal(distances,np.inf)
    reference=np.argsort(distances,axis=1)[:,:5]
    np.testing.assert_array_equal(idx,reference)
    # Reordering unique points preserves neighbour identities after alignment.
    perm=np.random.default_rng(20261109).permutation(len(X))
    permuted=nonself_neighbors(X[perm],5)
    np.testing.assert_array_equal(perm[permuted],idx[perm])
    return dict(duplicate_features_no_self_edges=True,independent_distance_reference=True,
                unique_distance_row_permutation=True,scope='Numerical graph invariants; ties beyond kth neighbour are not a global permutation certificate')
