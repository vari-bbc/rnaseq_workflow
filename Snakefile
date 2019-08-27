import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.4")


##### load config and sample sheets #####

configfile: "src/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"]).set_index("sample", drop=False)

contrasts = pd.read_table(config["contrasts"]).set_index("contrast", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

##### target rules #####

rule all:
    input:
        # symlink
        expand("raw_reads/{units.sample}_{units.unit}_{read}.fastq.gz", read=["R1","R2"], units=units.itertuples()),
        # # fastqc
        # expand("qc/fastqc/{units.sample}_{units.unit}_{read}_fastqc.html", read=["R1","R2"], units=units.itertuples()),
        # expand("qc/fastqc/{units.sample}_{units.unit}_{read}_fastqc.zip", read=["R1","R2"], units=units.itertuples()),
        # # Trim_Galore
        expand("trimmed_data/{units.sample}_{units.unit}_R1_val_1.fq.gz", units=units.itertuples()),
        expand("trimmed_data/{units.sample}_{units.unit}_R2_val_2.fq.gz", units=units.itertuples()),
        # # STAR alignment
        expand("analysis/star/{units.sample}_{units.unit}.Aligned.out.bam", units=units.itertuples()),
        expand("analysis/star/{units.sample}_{units.unit}.Log.out", units=units.itertuples()),
        # # sort_index_bam
        expand("analysis/star/{units.sample}_{units.unit}.sorted.bam", units=units.itertuples()),
        expand("analysis/star/{units.sample}_{units.unit}.sorted.bam.bai", units=units.itertuples()),
        # # multiqc
        "qc/multiqc_report.html",
        # # starMatrix
        "deliverables/UniquelyMappingRates.txt",
        "deliverables/UniquelyMappingReads.txt",
        "deliverables/starMatrix.txt",
        # counts
        "deliverables/counts.tsv",

        # # edger
        "deliverables/edgeR_shortReport.html",
        #"deliverables/edgeR_longReport.html",

##### load rules #####

include: "rules/common.smk"
include: "rules/symlink.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/count_matrix.smk"
include: "rules/qc.smk"
include: "rules/edger.smk"
