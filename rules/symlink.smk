# make list of sample names:
SAMPLES = samples["sample"].tolist()

# set up symlinks for files with names that don't fit standard convention
rule symlink:
    output:
        expand("raw_reads/{sample}_{unit}_R1.fastq.gz", sample=samples["sample"].tolist(), unit=set(units["unit"].tolist())),
        expand("raw_reads/{sample}_{unit}_R2.fastq.gz", sample=samples["sample"].tolist(), unit=set(units["unit"].tolist())),
    priority: 2
    shell:
        "Rscript src/symlinks.R"
