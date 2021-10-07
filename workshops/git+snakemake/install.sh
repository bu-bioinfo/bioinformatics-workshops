module load miniconda
conda create --name workshop seaborn scikit-learn statsmodels numba pytables
conda activate workshop
conda install -c conda-forge igraph python-igraph leidenalg imagemagick 
pip install scanpy snakemake scikit-misc louvain jinja2 pygraphviz pygments
