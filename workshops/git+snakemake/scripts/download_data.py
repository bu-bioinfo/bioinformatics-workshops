import scanpy as sc

def download_data(dataset):
    """Download and return specified dataset."""
    if dataset == 'pbmc3k':
        adata = sc.datasets.pbmc3k()
    elif dataset == 'paul':
        adata = sc.datasets.paul15()
    elif dataset == 'moignard':
        adata = sc.datasets.moignard15()
    else:
        raise ValueError(f"Unsupported dataset {dataset}.")
    return adata

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        adata = download_data(snakemake.params['dataset'])
        adata.write(snakemake.output['adata'])