# Bulk RNAseq Workflow

This workflow performs a differential expression analysis with STAR and edgeR - inspired by https://github.com/snakemake-workflows/rna-seq-star-deseq2

## Authors

* Dean Pettinga [@deanpettinga](https://github.com/deanpettinga)

## Usage

**NOTE** this workflow is optimized for HPC3 @ Van Andel Institute.

### Requirements

* You need to have an installation of [conda](https://docs.conda.io/en/latest/miniconda.html#linux-installers) in your `$PATH`
  * Snakemake is invoked from within a conda environment, so you must have your shell initialized for conda.
* You need to have [Singularity](https://sylabs.io/guides/3.4/user-guide/) installed and callable in your `$PATH`.
  * *Recommended*: by default, Singularity will cache containers in your `$HOME`. As user storage is limited on HPC3, please assign the environmental variables `$SINGULARITY_CACHEDIR` and `$SINGULARITY_TEMPDIR` to a different, non-limited directory to ensure that Singularity doesn't fail due to lack of cache space. E.g. (in your `~/.bash_profile`):
    ```
    # Set Singularity Cache to new location (default is $HOME/.singularity)
    export SINGULARITY_CACHEDIR=/path/to/singularity/cache/
    export SINGULARITY_TEMPDIR=/path/to/singularity/cache/
    ```

### Step 1: Installation

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a project directory on /secondary
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for your project.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).

### Step 2: Configure the workflow
* Modify the config, and any necessary sheets, e.g.:
  * src/units.tsv
    * **sample**        - ID of biological sample
    * **group**         - comparison group for DE contrast
    * **unit**          - description of sequencing unit ("lane1" or "flowcell1", etc.)
    * **fq1**           - name of read1 fastq
    * **fq2***           - name of read2 fastq
    * **strandedness**  - strandedness of library prep
      * typically reverse for Illumina
      * used to identify column to count reads from STAR output.
  * src/contrasts.tsv
    * **contrast**      - name for edgeR contrast
    * **groupRelative** - "group" defined in units.tsv
    * **groupBaseline** - other "group" defined in units.tsv
  * src/config.yaml
    * **PE_or_SE** - designates the sequencing protocol as paired-end or single-end. (e.g. "PE")
    * **ref**
      * **index** - absolute path to your STAR index directory
      * **annotation** - absolute path to your genome annotation.gtf
    * **species**
      * **short** - e.g. "hsapiens"
      * **long**  - e.g. "Homo sapiens"
    * **annotation** - [Bioconductor annotation package](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb) (e.g. "org.Hs.eg.db")
* Move your sequencing reads to `raw_reads/`

### Step 3: Test the workflow
Test your configuration by performing a dry-run via

    $ snakemake --use-conda -np

Execute from within your project directory as a PBS job using BBC nodes via

    $ qsub -q bbc src/run_snake.sh

This job script will produce DAG (.txt & .png) and .html with run stats for the workflow to be executed in `logs/runs/bulk_rnaseq-workflow_(TIME)`

### Step 4: Investigate Results
Review your results including the run stats:

* `logs/runs/rnaseq-workflow_(TIME).html`

and the differential expression results:
* `deliverables/edgeR_longReport.html`

### Step 5:
Now that you've successfully run your analysis, it's time to do some housekeeping.
* Commit any changes you've made to the repo and push the project-branch to your fork on github.
  * Optional: Merge back any valuable and generalizable changes to the [upstream repo] via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
  * Optional: Push results (plots/tables) to the remote branch on your fork.
  * Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
  * Optional: Delete the local clone/workdir to free space.
