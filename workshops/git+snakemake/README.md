# snakemake-workshop
This workshop should serve both as an introduction to `Snakemake` as well as `git` + `GitHub`. 

[readthedocs](https://programming-workshops.readthedocs.io/en/latest/workshops/04_git_cv/04_version_control_with_git_workshop.html)

## Workshop Overview
In this workshop you will use snakemake to create a workflow to download, process, and plot single-cell data.

## Getting Started
0. Confirm that you have read through the pre-workshop material and are familiar with the following concepts (if not, please do this now):
    - git status, git add, git commit, git pull, git push
    - what is a Snakefile
    - Open up the snakemake documentation and feel free to reference this as you work: <https://snakemake.readthedocs.io/en/stable/index.html>
1. Log on to the SCC
    - Hint: Start an interactive node
2. Fork the workshop repository to your own account: <https://github.com/BRITE-REU/snakemake-workshop>
3. Run each line of the [install script](install.sh) line-by-line on the command line
    - Check: Is the workshop environment activated? (try 'conda info')
4. Make the pipeline work
    - Hint: Navigate to the **Issues** board on your forked repository. The issues are numbered in the order in which you should work on them. When you finish an issue, mark it as closed. You can then view your progress in the **Projects** board.
    - Hint: Use git to commit and push your changes as you work

## Some Hints
- To test your work as you go, we suggest you comment out all rules that follow the current rule you are implementing. Then when you call the Snakfile, specify which output file you want to create. 
