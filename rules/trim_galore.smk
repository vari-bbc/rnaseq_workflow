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
        "trimmed_data/{sample}-R1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}-R2_val_2.fq.gz",
        "trimmed_data/{sample}-R2.fastq.gz_trimming_report.txt"
    params:
        extra = "-q 20"
    log:
        "logs/trim/trim_galore.{sample}.log"
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
    wrapper:
        #"0.31.1/bio/trim_galore/pe"
        "file:wrappers/trim_galore_se"
