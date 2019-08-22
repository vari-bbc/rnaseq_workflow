# Bulk RNAseq Workflow

This workflow performs a differential expression analysis with STAR and edgeR - strongly inspired by https://github.com/snakemake-workflows/rna-seq-star-deseq2

## Authors

* Dean Pettinga (@deanpettinga), https://github.com/deanpettinga

## Usage

**NOTE** this workflow is optimized for HPC3 @ Van Andel Institute.

### Step 1: Installation

make sure you are running [![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.4-green.svg)](https://snakemake.bitbucket.io)

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a project directory on /secondary
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for your project.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).

### Step 2: Configure the workflow
* Modify the config, and any necessary sheets, e.g.:
  * src/samples.tsv
  * src/contrasts.tsv
  * src/config.yaml
  * src/cluster.json
    * **NOTE** - default reference is ensembl hg38 gencode annotation
* Move your sequencing reads to `raw_reads/`

### Step 3: Test the workflow
Test your configuration by performing a dry-run via

    snakemake --use-conda -np

Execute as from within your project directory as a PBS job using BBC nodes via

    qsub -q bbc /src/run_snake.sh

This job script will produce DAG (.txt & .png) and .html with run stats for the workflow to be executed in `runs/bulk_rnaseq_workflow_(TIME)`

### Step 4: Investigate Results
Review your results including the run stats:

* `runs/bulk_rnaseq_workflow_(TIME).html`

and the differential expression results:
* `deliverables/edgeR_longReport.html`
* `deliverables/edgeR_shortReport.html`

### Step 5:
Now that you've successfully run your analysis, it's time to do some housekeeping.
* Commit any changes you've made to the repo and push the project-branch to your fork on github.
  * Optional: Merge back any valuable and generalizable changes to the [upstream repo] via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
  * Optional: Push results (plots/tables) to the remote branch on your fork.
  * Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
  * Optional: Delete the local clone/workdir to free space.
