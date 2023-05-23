# Bulk RNAseq Workflow

This workflow read trimming and read alignment with STAR, producing a read count smatrix that can be used for differential expression analysis.

* [Bulk RNAseq Workflow](#bulk-rnaseq-workflow)
   * [Usage](#usage)
      * [Step 1: Configure the workflow](#step-1-configure-the-workflow)
      * [Step 2: Test and run the workflow](#step-2-test-and-run-the-workflow)
      * [Optional: Variant calling](#optional-variant-calling)


## Usage

**NOTE** this workflow is optimized for HPC3 @ Van Andel Institute.


### Step 1: Configure the workflow
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
  * bin/config.yaml
    * **PE_or_SE** - designates the sequencing protocol as paired-end or single-end. (e.g. "PE")
    * **ref**
      * **index** - absolute path to your STAR index directory
      * **annotation** - absolute path to your genome annotation.gtf
* Move your sequencing reads to `raw_data/`

### Step 2: Test and run the workflow
Test your configuration by performing a dry-run via

    $ snakemake -np

Execute from within your project directory as a SLURM job.

    $ sbatch bin/run_snake.sh

This job script will produce DAG (.txt & .png) and .html with run stats for the workflow to be executed in `logs/runs/bulk_rnaseq-workflow_(TIME)`


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
