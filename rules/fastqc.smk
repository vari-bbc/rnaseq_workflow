rule fastqc_PE:
    input:
        "raw_reads/{sample}_{unit}_R{read}.fastq.gz",
    output:
        html="qc/fastqc/{sample}_{unit}_R{read}_fastqc.html",
        zip="qc/fastqc/{sample}_{unit}_R{read}_fastqc.zip",
    params: ""
    wrapper:
        "0.35.1/bio/fastqc"

rule fastqc_SE:
    input:
        "raw_reads/{sample}-{unit}.fastq.gz",
    output:
        html="qc/fastqc/{sample}-{unit}_fastqc.html",
        zip="qc/fastqc/{sample}-{unit}_fastqc.zip",
    params: ""
    wrapper:
        "0.35.1/bio/fastqc"
