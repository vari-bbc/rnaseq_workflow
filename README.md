# Bulk RNAseq Workflow


* [Bulk RNAseq Workflow](#bulk-rnaseq-workflow)
   * [Usage](#usage)
      * [Step 1: Configure the workflow](#step-1-configure-the-workflow)
      * [Step 2: Test and run the workflow](#step-2-test-and-run-the-workflow)

## Usage

**NOTE** this workflow is optimized for the HPC @ Van Andel Institute.


### Step 1: Configure the workflow
* Move your sequencing reads to `raw_data/`

* Modify the config and samplesheet:
  * config/samplesheet/units.tsv; To make a template based in the files in `raw_data/`, run `./make_units_template.sh`.
    * **sample**        - ID of biological sample; Must be unique.
    * **group**         - Experimental group 
    * **fq1**           - name of read1 fastq
    * **fq2***          - name of read2 fastq

  * config/config.yaml

### Step 2: Test and run the workflow
Test your configuration by performing a dry-run via

    $ snakemake -np

Execute from within your project directory as a SLURM job.

    $ sbatch bin/run_snake.sh



