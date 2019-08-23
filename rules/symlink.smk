# make list of sample names:
SAMPLES = samples["sample"].tolist()

# set up symlinks for files with names that don't fit standard convention
rule symlink:
    output:
        expand("raw_reads/{sample}_1.fastq.gz", sample=SAMPLES),
        expand("raw_reads/{sample}_2.fastq.gz",  sample=SAMPLES)
    priority: 2
    shell:
        "Rscript src/symlinks.R"
