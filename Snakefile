import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.4")


##### load config and sample sheets #####

configfile: "src/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

units = pd.read_table(config["units"]).set_index("sample", drop=False)
validate(units, schema="schemas/units.schema.yaml")

contrasts = pd.read_table(config["contrasts"]).set_index("contrast", drop=False)
validate(contrasts, schema="schemas/contrasts.schema.yaml")


##### target rules #####

rule all:
    input:
        # symlink
        # expand("raw_reads/{units.sample}_{units.unit}_R1.fastq.gz", units=units.itertuples()),
        # expand("raw_reads/{units.sample}_{units.unit}_R2.fastq.gz", units=units.itertuples()),
        # fastqc
        # expand("qc/fastqc/{units.sample}_{units.unit}_R1_fastqc.html", units=units.itertuples()),
        # expand("qc/fastqc/{units.sample}_{units.unit}_R2_fastqc.html", units=units.itertuples()),
        # Trim_Galore
        # expand("trimmed_data/{units.sample}_{units.unit}_R1_val_1.fq.gz", units=units.itertuples()),
        # expand("trimmed_data/{units.sample}_{units.unit}_R2_val_2.fq.gz", units=units.itertuples()),
        # STAR alignment
        # expand("analysis/star/{units.sample}_{units.unit}.Aligned.out.bam", units=units.itertuples()),
        # expand("analysis/star/{units.sample}_{units.unit}.Log.out", units=units.itertuples()),
        # sort_index_bam
        # expand("analysis/star/{units.sample}_{units.unit}.sorted.bam", units=units.itertuples()),
        # expand("analysis/star/{units.sample}_{units.unit}.sorted.bam.bai", units=units.itertuples()),
        # count_matrix
        # "deliverables/counts.tsv",
        # "deliverables/UniquelyMappingRates.txt",
        # "deliverables/UniquelyMappingReads.txt",
        # "deliverables/starMatrix.txt",
        # multiqc
        # "qc/multiqc_report.html",
        # edger
        "edgeR_longReport.html",

##### load rules #####

include: "rules/common.smk"
include: "rules/symlink.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/count_matrix.smk"
include: "rules/qc.smk"
include: "rules/edger.smk"
