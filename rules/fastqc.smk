rule fastqc_PE:
    input:
        "raw_reads/{sample}-R{read}.fastq.gz",
    output:
        html="qc/fastqc/{sample}-R{read}_fastqc.html",
        zip="qc/fastqc/{sample}-R{read}_fastqc.zip",
    params: ""
    wrapper:
        "0.35.1/bio/fastqc"

rule fastqc_SE:
    input:
        "raw_reads/{sample}-SE.fastq.gz",
    output:
        html="qc/fastqc/{sample}-SE_fastqc.html",
        zip="qc/fastqc/{sample}-SE_fastqc.zip",
    params: ""
    wrapper:
        "0.35.1/bio/fastqc"
