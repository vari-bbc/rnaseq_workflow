import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.4.4")

##### load config and sample sheets #####

configfile: "src/config.yaml"
validate(config, schema="schemas/config.schema.yaml")

units = pd.read_table(config["units"]).set_index("sample", drop=False)
validate(units, schema="schemas/units.schema.yaml")

contrasts = pd.read_table(config["contrasts"]).set_index("name", drop=False)
validate(contrasts, schema="schemas/contrasts.schema.yaml")

##### target rules #####

rule all:
    input:
        # mergeLanesAndRename
            #PE
        # expand("raw_reads/{sample}-R1.fastq.gz",sample=set(units["sample"].tolist())),
        # expand("raw_reads/{sample}-R2.fastq.gz",sample=set(units["sample"].tolist())),
            #SE
        #expand("raw_reads/{units.sample}.fastq.gz", units=units.itertuples()),
        # fastqc
        # expand("qc/fastqc/{units.sample}_fastqc.html", units=units.itertuples()),
        # expand("qc/fastqc/{units.sample}-R1_fastqc.html", units=units.itertuples()),
        # expand("qc/fastqc/{units.sample}-R2_fastqc.html", units=units.itertuples()),
        # Trim_Galore
        # expand("trimmed_data/{units.sample}-{units.unit}_trimmed.fq.gz", units=units.itertuples()), #SE
        # expand("trimmed_data/{units.sample}_{units.unit}_R1_val_1.fq.gz", units=units.itertuples()),
        # expand("trimmed_data/{units.sample}_{units.unit}_R2_val_2.fq.gz", units=units.itertuples()),
        # STAR alignment
        #expand("analysis/star/{units.sample}.Aligned.sortedByDoord.out.bam", units=units.itertuples()),
        #expand("analysis/star/{units.sample}.Log.out", units=units.itertuples()),
        # te_count
        # expand("analysis/TEcount/{units.sample}.cntTable", units=units.itertuples()),
        # multiQC
        "qc/multiqc_report.html",
        #"rules/diffExp.html"
        # mergeCounts
        # "deliverables/GeneTEcounts.tsv",
        "deliverables/GeneCounts.tsv"

##### load rules #####
include: "rules/mergeLanesAndRename.smk"
include: "rules/fastqc.smk"
include: "rules/trim_galore.smk"
include: "rules/makeIndex.smk"
include: "rules/STAR.smk"
include: "rules/starTEs.smk"
include: "rules/TEcount.smk"
include: "rules/mergeTEcounts.smk"
include: "rules/multiqc.smk"
include: "rules/edger.smk"
