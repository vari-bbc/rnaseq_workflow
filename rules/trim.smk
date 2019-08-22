def get_fastq(wildcards):
    return expand("raw_reads/"+samples.loc[wildcards.sample, ["sample"]]+"_{read}.fastq.gz", read=[1,2])

rule trim_galore_pe:
    input:
        # get_fastq
        "raw_reads/{sample}_1.fastq.gz",
        "raw_reads/{sample}_2.fastq.gz"
    output:
        "trimmed_data/{sample}_1_val_1.fq.gz",
        "trimmed_data/{sample}_1.fastq.gz_trimming_report.txt",
        "trimmed_data/{sample}_2_val_2.fq.gz",
        "trimmed_data/{sample}_2.fastq.gz_trimming_report.txt"
    params:
        extra = "--retain_unpaired -q 20"
    log:
        "logs/trim_galore.{sample}.log"
    conda:
        "../envs/trimgalore.yaml"
    wrapper:
        "0.31.1/bio/trim_galore/pe"
