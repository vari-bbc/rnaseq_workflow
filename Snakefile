import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.10.0")

##### load config and sample sheets #####

configfile: "src/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

units = pd.read_table(config["units"]).set_index("sample", drop=False)
#validate(units, schema="schemas/units.schema.yaml")

contrasts = pd.read_table(config["contrasts"]).set_index("name", drop=False)
#validate(contrasts, schema="schemas/contrasts.schema.yaml")

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
        # expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        # expand("analysis/star/{units.sample}.Log.out", units=units.itertuples()),
        # multiQC
        #"qc/multiqc_report.html",
        "src/diffExp.html"

rule mergeLanesAndRename_SE:
    input:
    output:
        "raw_reads/{sample}-SE.fastq.gz"
    log:
        "logs/mergeLanesAndRename/mergeLanesAndRename_SE-{sample}.log"
    threads: 1
    resources:
        mem_gb = 16
    conda:
        "envs/R.yaml"
    script:
        "src/mergeLanesAndRename.R"

rule mergeLanesAndRename_PE:
    input:
    output:
        "raw_reads/{sample}-R1.fastq.gz",
        "raw_reads/{sample}-R2.fastq.gz"
    log:
        "logs/mergeLanesAndRename/mergeLanesAndRename_PE-{sample}.log"
    threads: 1
    resources:
        mem_gb = 16
    conda:
        "envs/R.yaml"
    script:
        "src/mergeLanesAndRename.R"

rule fastqc_PE:
    input:
        "raw_reads/{sample}-R{read}.fastq.gz",
    output:
        html="qc/fastqc/{sample}-R{read}_fastqc.html",
        zip="qc/fastqc/{sample}-R{read}_fastqc.zip",
    params: ""
    threads: 1
    resources:
        mem_gb = 64
    wrapper:
        "file:wrappers/fastqc"

rule fastqc_SE:
    input:
        "raw_reads/{sample}-SE.fastq.gz",
    output:
        html="qc/fastqc/{sample}-SE_fastqc.html",
        zip="qc/fastqc/{sample}-SE_fastqc.zip",
    params: ""
    threads: 1
    resources:
        mem_gb = 64
    wrapper:
        "file:wrappers/fastqc"

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "raw_reads/{sample}-SE.fastq.gz".format(**wildcards)
        fastqc_html = "qc/fastqc/{sample}-SE_fastqc.html".format(**wildcards)
        fastqc_zip = "qc/fastqc/{sample}-SE_fastqc.zip".format(**wildcards)
        return [reads,fastqc_html,fastqc_zip]
    elif config["PE_or_SE"] == "PE":
        R1 = "raw_reads/{sample}-R1.fastq.gz".format(**wildcards)
        R2 = "raw_reads/{sample}-R2.fastq.gz".format(**wildcards)
        fastqc_html_1 = "qc/fastqc/{sample}-R1_fastqc.html".format(**wildcards)
        fastqc_html_2 = "qc/fastqc/{sample}-R2_fastqc.html".format(**wildcards)
        fastqc_zip_1 = "qc/fastqc/{sample}-R1_fastqc.zip".format(**wildcards)
        fastqc_zip_2 = "qc/fastqc/{sample}-R2_fastqc.zip".format(**wildcards)
        return [R1,R2,fastqc_html_1,fastqc_html_2,fastqc_zip_1,fastqc_zip_2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        "trimmed_data/{sample}-R1_val_1.fq.gz",
        "trimmed_data/{sample}-R1_val_1_fastqc.html",
        "trimmed_data/{sample}-R1_val_1_fastqc.zip",
        "trimmed_data/{sample}-R1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}-R2_val_2.fq.gz",
        "trimmed_data/{sample}-R2_val_2_fastqc.html",
        "trimmed_data/{sample}-R2_val_2_fastqc.zip",
        "trimmed_data/{sample}-R2.fastq.gz_trimming_report.txt"
    params:
        extra = "-q 20"
    log:
        "logs/trim/trim_galore.{sample}.log"
    threads: 4
    resources:
        mem_gb = 64
    wrapper:
        #"0.31.1/bio/trim_galore/pe"
        "file:wrappers/trim_galore_pe"

rule trim_galore_SE:
    input:
        trim_galore_input
    output:
        "trimmed_data/{sample}-SE_trimmed.fq.gz",
        "trimmed_data/{sample}-SE.fastq.gz_trimming_report.txt",
    params:
        extra = "--illumina -q 20"
    log:
        "logs/trim/trim_galore.{sample}.log"
    threads: 4
    resources:
        mem_gb = 64
    wrapper:
        #"0.31.1/bio/trim_galore/pe"
        "file:wrappers/trim_galore_se"

def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="trimmed_data/{sample}-SE_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "trimmed_data/{sample}-R1_val_1.fq.gz".format(**wildcards)
        fq2 = "trimmed_data/{sample}-R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule STAR:
    input:
        STAR_input
    output:
        # see STAR manual for additional output files
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
        "analysis/star/{sample}.Log.final.out",
        "analysis/star/{sample}.Log.out",
        "analysis/star/{sample}.Log.progress.out",
        "analysis/star/{sample}.ReadsPerGene.out.tab",
        "analysis/star/{sample}.SJ.out.tab",
        directory("analysis/star/{sample}._STARgenome"),
        directory("analysis/star/{sample}._STARpass1"),

    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts "
    threads: 8
    resources:
        mem_gb = 120
    wrapper:
        "file:wrappers/star"

multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("qc/fastqc/{units.sample}-SE_fastqc.html", units=units.itertuples()))
    multiqc_input.append(expand("qc/fastqc/{units.sample}-SE_fastqc.zip", units=units.itertuples()))
    multiqc_input.append(expand("trimmed_data/{units.sample}-SE.fastq.gz_trimming_report.txt", units=units.itertuples()))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("qc/fastqc/{units.sample}-{read}_fastqc.html", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("qc/fastqc/{units.sample}-{read}_fastqc.zip", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("trimmed_data/{units.sample}-{read}.fastq.gz_trimming_report.txt", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.html", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.zip", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))

rule multiqc:
    input:
        multiqc_input
    params:
        # skip the pass1 from STAR
        "--ignore '*._STARpass1/*'"
    output:
        "qc/multiqc_report.html",
        "qc/multiqc_report_data/multiqc.log",
        "qc/multiqc_report_data/multiqc_cutadapt.txt",
        "qc/multiqc_report_data/multiqc_fastqc.txt",
        "qc/multiqc_report_data/multiqc_general_stats.txt",
        "qc/multiqc_report_data/multiqc_sources.txt",
        "qc/multiqc_report_data/multiqc_star.txt",
    log:
        "logs/multiqc.log"
    threads: 1
    resources:
        mem_gb = 16
    wrapper:
        "file:wrappers/multiqc"

rule edgeR:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        "qc/multiqc_report.html" # require multiQC to be run before this analysis
    output:
        "src/diffExp.html"
    log:
        "logs/edgeR.log"
    singularity:
        "docker://dpettinga/bbcrna:latest"
    threads: 1
    resources:
        mem_gb = 16
    script:
        "src/diffExp.Rmd"
