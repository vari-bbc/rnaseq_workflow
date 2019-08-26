rule trim_galore_pe:
    input:
        # get_fastq
        "raw_reads/{sample}_{unit}_R1.fastq.gz",
        "raw_reads/{sample}_{unit}_R2.fastq.gz"
    output:
        "trimmed_data/{sample}_{unit}_1_val_1.fq.gz",
        "trimmed_data/{sample}_{unit}_1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}_{unit}_2_val_2.fq.gz",
        "trimmed_data/{sample}_{unit}_2.fastq.gz_trimming_report.txt"
    params:
        extra = "--retain_unpaired -q 20"
    log:
        "logs/trim_galore.{sample}_{unit}.log"
    conda:
        "../envs/trimgalore.yaml"
    wrapper:
        "0.31.1/bio/trim_galore/pe"
