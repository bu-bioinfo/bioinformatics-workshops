# snakemake-workshop
Worklflow managers are good mkay

## Workshop Overview
In this workshop you will use snakemake to create a workflow to download, process, and plot single-cell data. (sea urchins)

## Getting Started
1. Log on to the SCC
    - Hint: Start an interactive node
2. Fork the workshop repository to your own account: <https://github.com/BRITE-REU/snakemake-workshop>
3. Run each line of the [install script](install.sh) line-by-line on the command line
    - Check: Is the workshop environment activated? (try 'conda info')
4. Make the pipeline work
    - Hint: Navigate to the **Issues** board on your forked repository. The issues are numbered in the order in which you should work on them. When you finish an issue, mark it as closed. You can then view your progress in the **Projects** board.
    - Hint: Use git to commit and push your changes as you work

TODO:

1. Finish Pipeline (no R)
2. Test Pipeline
3. Generalize / finalize pipeline for multiple datasets
5. Remove everything expect rule names + descriptions from Snakemake
6. Remove all snakemake handling in Python + R Scripts
7. Create "TODO" list for workshop using issues in github
  - Link to workshop project for automated khan ban board
  - Start issues as only implementing rule for single dataset (pbmc3k)
  - Introduce config file to hold parameters
  - Generalize rules to handle multipe wildcards
  - Generate report 
8. Merge branch to master. 
  - Previous commit now has solutions
  - Current commit now has issues detaliing what students should do, they can fork ezpz
