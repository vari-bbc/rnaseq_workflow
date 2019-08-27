rule trim_galore_pe:
    input:
        # this ensures that fastqc occurs first
        "qc/fastqc/{sample}_{unit}_R1_fastqc.html",
        "qc/fastqc/{sample}_{unit}_R1_fastqc.zip"
        "qc/fastqc/{sample}_{unit}_R2_fastqc.html",
        "qc/fastqc/{sample}_{unit}_R2_fastqc.zip"
        # get_fastq
        "raw_reads/{sample}_{unit}_R1.fastq.gz",
        "raw_reads/{sample}_{unit}_R2.fastq.gz"
    output:
        "trimmed_data/{sample}_{unit}_R1_val_1.fq.gz",
        "trimmed_data/{sample}_{unit}_R1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}_{unit}_R2_val_2.fq.gz",
        "trimmed_data/{sample}_{unit}_R2.fastq.gz_trimming_report.txt"
    params:
        extra = "--retain_unpaired -q 20"
    log:
        "logs/trim/trim_galore.{sample}_{unit}.log"
    conda:
        "../envs/trimgalore.yaml"
    wrapper:
        "0.31.1/bio/trim_galore/pe"
