# Pipeline to download, preprocess, cluster, and plot cells in example scRNAseq data
# @authors: Dakota Hawkins and Emma Briars


# This should call the download_data.py Python script to download data.
# You will need to tell the script which dataset you want to download, and the
# file name to save the dataset under. The expected extension for the output file
# is an '.h5ad' file. In snakemake, it is good practice to write data files
# produced by each rule to unique directories.
# allowable: paul, moignard, pbmc3k
rule download_data:

# This rule should preprocess downloaded data by calling the `preprocess.py`
# Python script.
# The rule should read in the raw dataset, and filter cells + genes based off
# of provided parameters as explained in the github issue. The rule should write
# the newly processed data to a new `.h5ad` file. 
rule preprocess_data:


# This rule should cluster cells using the `cluster_cells.py` script.
# rule should read in preprocessed data, and clsuter cells according to user-
# provided parameters, k and resolution. The rule should produce 3 output files
# as csvs: a count matrix, a cell metadata table, and a gene metadata table.
rule cluster_cells:


# This rule should plot clusters on UMAP projections using ggplot and the 
# plot_cells.R script. 
# the rule should read in the previously generated csvs, color cells according
# to values in a user-specified column, and create a .png file containing the plot
rule plot_clusters:
