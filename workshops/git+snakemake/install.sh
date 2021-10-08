module load miniconda
conda create --name workshop seaborn scikit-learn statsmodels numba pytables pip
conda activate workshop
conda install -c conda-forge igraph python-igraph leidenalg imagemagick 
pip install "scanpy[louvain]" snakemake scikit-misc jinja2 pygraphviz pygments
