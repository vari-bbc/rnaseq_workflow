# def get_mergeLanesAndRename_input(wildcards):
#     if config["PE_or_SE"]=="SE":
#         return expand("raw_reads/{fq1}",**wildcards)
#     elif config["PE_or_SE"]=="PE":
#         read_list = [expand("raw_reads/{sample}-R1.fastq.gz",sample=set(units["sample"].tolist())),
#         expand("raw_reads/{sample}-R2.fastq.gz",sample=set(units["sample"].tolist()))]
#         flat_list = [read for read_set in read_list for read in read_set]
#         return flat_list

rule mergeLanesAndRename_SE:
    input:
    output:
        "raw_reads/{sample}-SE.fastq.gz"
    log:
        "logs/mergeLanesAndRename_SE.{sample}.log"
    conda:
        "../envs/R.yaml"
    log:
        "logs/mergeLanesAndRename_SE-{sample}.log"
    script:
        "mergeLanesAndRename.R"

rule mergeLanesAndRename_PE:
    input:
    output:
        "raw_reads/{sample}-R1.fastq.gz",
        "raw_reads/{sample}-R2.fastq.gz"
    log:
        "logs/mergeLanesAndRename_PE-{sample}.log"
    conda:
        "../envs/R.yaml"
    script:
        "mergeLanesAndRename.R"
