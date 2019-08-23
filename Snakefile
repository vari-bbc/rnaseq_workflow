import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.4")


##### load config and sample sheets #####

configfile: "src/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

contrasts = pd.read_table(config["contrasts"]).set_index("contrast", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

# make list of sample names:
SAMPLES = samples["sample"].tolist()
##### target rules #####

rule all:
    input:
        # symlink
        expand("raw_reads/{samples.sample}_{read}.fastq.gz", read=["1","2"], samples=samples.itertuples()),
        # fastqc_raw
        expand("qc/fastqc/raw_reads/{samples.sample}_{read}_fastqc.html", read=["1","2"], samples=samples.itertuples()),
        expand("qc/fastqc/raw_reads/{samples.sample}_{read}_fastqc.zip", read=["1","2"], samples=samples.itertuples()),
        # Trim_Galore
        expand("trimmed_data/{samples.sample}_1_val_1.fq.gz", samples=samples.itertuples()),
        expand("trimmed_data/{samples.sample}_2_val_2.fq.gz", samples=samples.itertuples()),
        # STAR alignment
        expand("star/{samples.sample}.Aligned.out.bam", samples=samples.itertuples()),
        expand("star/{samples.sample}.Log.out", samples=samples.itertuples()),
        # sort_index_bam
        expand("star/{samples.sample}.sorted.bam", samples=samples.itertuples()),
        expand("star/{samples.sample}.sorted.bam.bai", samples=samples.itertuples()),
        # multiqc
        "qc/multiqc_report.html",
        # starMatrix
        "deliverables/UniquelyMappingRates.txt",
        "deliverables/UniquelyMappingReads.txt",
        "deliverables/starMatrix.txt",
        # edger
        #"deliverables/edgeR_shortReport.html",
        "deliverables/edgeR_longReport.html",

##### load rules #####

include: "rules/common.smk"
include: "rules/symlink.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/count_matrix.smk"
include: "rules/qc.smk"
include: "rules/edger.smk"
