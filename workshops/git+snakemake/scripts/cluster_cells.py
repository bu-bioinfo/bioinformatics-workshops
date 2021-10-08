import numpy as np
import pandas as pd
import scanpy as sc

def cluster_cells(adata, k, res):
    """
    Cluster cells using louvain community detection

    params
    ------
        adata : sc.AnnData
            dataset to cluster
        k : int
            number of neighbors to use when constructing the neighbor graph.
        res : float
            resolution parameter for louvain.
    returns
    -------
    sc.AnnData
        Clustered dataset with UMAP projections.
    """
    sc.pp.pca(adata, n_comps=min(adata.shape[1] - 1, 50))
    sc.pp.neighbors(adata, n_neighbors=k)
    sc.tl.louvain(adata, resolution=res)
    sc.tl.umap(adata)
    adata.obs['umap1'] = adata.obsm['X_umap'][:, 0]
    adata.obs['umap2'] = adata.obsm['X_umap'][:, 1]
    return adata

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        # read in data using sc.read()
        adata = sc.read()
        out = cluster_cells(adata)
        # writes count matrix to a csv file
        np.savetxt(, out.X, delimiter=',')
        # writes cell observation/metadata to csv file
        out.obs.to_csv()
        # writes gene metadata to csv file
        out.var.to_csv()
