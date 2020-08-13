# Bulk RNAseq Workflow

This workflow performs a differential expression analysis with STAR and edgeR

* [Bulk RNAseq Workflow](#bulk-rnaseq-workflow)
   * [Authors](#authors)
   * [Usage](#usage)
      * [Step 1: Installation](#step-1-installation)
      * [Step 2: Configure the workflow](#step-2-configure-the-workflow)
      * [Step 3: Test the workflow](#step-3-test-the-workflow)
      * [Step 4: Investigate Results](#step-4-investigate-results)
      * [Step 5: 'Housekeeping'](#step-5-housekeeping)
      * [Optional: Variant calling](#optional-variant-calling)

## Authors

* Dean Pettinga [@deanpettinga](https://github.com/deanpettinga)


## Usage

**NOTE** this workflow is optimized for HPC3 @ Van Andel Institute.

### Step 1: Installation

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a project directory on /secondary
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for your project.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).

### Step 2: Configure the workflow
* Modify the config and any necessary sheets, e.g.:
  * bin/units.tsv
    * **sample**        - ID of biological sample
    * **group**         - comparison group for DE contrast
    * **unit**          - description of sequencing unit ("lane1" or "flowcell1", etc.)
    * **fq1**           - name of read1 fastq
    * **fq2***           - name of read2 fastq
    * **strandedness**  - strandedness of library prep
      * typically reverse for Illumina
      * used to identify column to count reads from STAR output.
  * bin/contrasts.tsv
    * **contrast**      - name for edgeR contrast
    * **groupRelative** - "group" defined in units.tsv
    * **groupBaseline** - other "group" defined in units.tsv
  * bin/config.yaml
    * **PE_or_SE** - designates the sequencing protocol as paired-end or single-end. (e.g. "PE")
    * **ref**
      * **index** - absolute path to your STAR index directory
      * **annotation** - absolute path to your genome annotation.gtf
    * **species**
      * **short** - e.g. "hsapiens"
      * **long**  - e.g. "Homo sapiens"
    * **annotation** - [Bioconductor annotation package](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb) (e.g. "org.Hs.eg.db")
* Move your sequencing reads to `raw_data/`

### Step 3: Test the workflow
Test your configuration by performing a dry-run via

    $ snakemake --use-conda -np

Execute from within your project directory as a PBS job using BBC nodes via

    $ qsub -q bbc bin/run_snake.sh

This job script will produce DAG (.txt & .png) and .html with run stats for the workflow to be executed in `logs/runs/bulk_rnaseq-workflow_(TIME)`

### Step 4: Investigate Results
Review your results including the run stats:

* `logs/runs/rnaseq-workflow_(TIME).html`

and the differential expression results:
* `deliverables/edgeR_longReport.html`

### Step 5: 'Housekeeping'
Now that you've successfully run your analysis, it's time to do some housekeeping.
* Commit any changes you've made to the repo and push the project-branch to your fork on github.
  * Optional: Merge back any valuable and generalizable changes to the [upstream repo] via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
  * Optional: Push results (plots/tables) to the remote branch on your fork.
  * Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
  * Optional: Delete the local clone/workdir to free space.


### Optional: Variant calling
Besides running the differential expression portion of the pipeline, it is possible to also do variant calling based on the GATK "Best Practices". Importantly, the read groups are added post-alignment to the STAR BAM files created in the DE part of the pipeline, so align files separately if you want to assign distinct reads groups.

1. Create a tab-separated file called `bin/variant_calling_units.tsv`. These may consist of all or a subset of the files used for DE analysis. There are four columns:

    i. unit: This corresponds the 'sample' prefix in `{sample}.Aligned.sortedByCoord.out.bam` from the STAR alignments and is used to code the `RG:ID` tag.

    ii. sample: This is used to merge same-sample 'units' together at the MarkDuplicates step. Specifically, this indicates the value for the `RG:SM` tag.

    iii. library: Indicates the value for the `RG:LB` tag. Important for the MarkDuplicates step.

    iv. platform_unit: Indicates the value for the `RG:PU` tag. According to GATK docs, this should be {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. This is important for the BQSR step; if fastq files were merged in the DE part of the pipeline, then you cannot correct for any potential batch effects due to lane, flow cell etc.

2. In the `config.yaml` file:

    i. Set `call_variants` to `true`.

    ii. Ensure that valid paths are provided for `known_snps`, `known_indels` and `fai`. These should be two vcf/vcf.gz files and a .fai file. For hg38, use snpDB and "known_indels" from the GATK resource bundle.

3. For this to work, you must specify the .fai file first in step 2. You must also ensure that your default python has access to the snakemake library. If running on the VAI HPC, you may use the python installed using venv for the environment module of snakemake. Run `python ./group_chroms.py`. This will create `grouped_contigs.tsv` which is used to run GATK in parallel in a per-contig fashion.

4. Do a dry run to check if the rules are being correctly triggered. Run the pipeline as usual.
