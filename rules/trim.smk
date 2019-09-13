def get_fastq(wildcards):
    x = units.loc[lambda df: df['sample'] == wildcards.sample].loc[lambda df: df['unit'] == wildcards.unit,("fq1","fq2")].values.tolist()[0]
    x = [read for read in x if str(read) != 'nan']
    return x

rule trim_galore_PE:
    input:
        reads=["raw_reads/{sample}_{unit}_R1.fastq.gz","raw_reads/{sample}_{unit}_R2.fastq.gz"],
        # this ensures that fastqc occurs before this step
        fastqc_html=["qc/fastqc/{sample}_{unit}_R1_fastqc.html","qc/fastqc/{sample}_{unit}_R2_fastqc.html"],
        fastqc_zip=["qc/fastqc/{sample}_{unit}_R1_fastqc.zip","qc/fastqc/{sample}_{unit}_R2_fastqc.zip"],
    output:
        "trimmed_data/{sample}_{unit}_R1_val_1.fq.gz",
        "trimmed_data/{sample}_{unit}_R1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}_{unit}_R2_val_2.fq.gz",
        "trimmed_data/{sample}_{unit}_R2.fastq.gz_trimming_report.txt"
    params:
        extra = "-q 20"
    log:
        "logs/trim/trim_galore.{sample}_{unit}.log"
    wrapper:
        #"0.31.1/bio/trim_galore/pe"
        "file:wrappers/trim_galore_pe"

rule trim_galore_SE:
    input:
        reads="raw_reads/{sample}-{unit}.fastq.gz",
        # this ensures that fastqc occurs before this step
        fastqc_html="qc/fastqc/{sample}-{unit}_fastqc.html",
        fastqc_zip="qc/fastqc/{sample}-{unit}_fastqc.zip",
    output:
        "trimmed_data/{sample}-{unit}_trimmed.fq.gz",
        "trimmed_data/{sample}-{unit}.fastq.gz_trimming_report.txt",
    params:
        extra = "--illumina -q 20"
    log:
        "logs/trim/trim_galore.{sample}_{unit}.log"
    wrapper:
        #"0.31.1/bio/trim_galore/pe"
        "file:wrappers/trim_galore_se"
