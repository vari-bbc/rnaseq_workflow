# set up symlinks for files with names that don't fit standard convention
rule symlink:
    output:
        expand("raw_reads/{sample}_{unit}_R1.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
        expand("raw_reads/{sample}_{unit}_R2.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
    priority: 2
    shell:
        "Rscript src/symlink.R"
