rule mergeLanesAndRename_SE:
    input:
        expand("raw_reads/{fq1}", fq1=units["fq1"]),
    output:
        expand("raw_reads/{sample}.fastq.gz",sample=set(units["sample"].tolist())),
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript src/mergeLanesAndRename.R"

rule mergeLanesAndRename_PE:
    input:
        expand("raw_reads/{fq1}", fq1=units["fq1"]),
        expand("raw_reads/{fq2}", fq2=units["fq2"]),
    output:
        expand("raw_reads/{sample}-R1.fastq.gz",sample=set(units["sample"].tolist())),
        expand("raw_reads/{sample}-R2.fastq.gz",sample=set(units["sample"].tolist())),
    conda:
        "../envs/R.yaml"
    shell:
        "Rscript src/mergeLanesAndRename.R"
