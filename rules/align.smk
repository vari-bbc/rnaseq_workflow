def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("trimmed_data/{sample}_{read_trimmed}.fq.gz",
                          read_trimmed=["1_val_1", "2_val_2"], **wildcards)
        # single end sample
        return "trimmed_data/{sample}.fq.gz".format(**wildcards)


rule align:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = "trimmed_data/{sample}_1_val_1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "trimmed_data/{sample}_2_val_2.fq.gz"
    output:
        # see STAR manual for additional output files
        "star/{sample}.ReadsPerGene.out.tab",
        "star/{sample}.Aligned.out.bam",
        "star/{sample}.Log.out",
        "star/{sample}.Log.progress.out",
        "star/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts"
    threads: 24
    wrapper:
        "file:wrappers/star"

rule sort_index_bam:
    input:
        "star/{sample}.Aligned.out.bam",
    output:
        sorted = "star/{sample}.sorted.bam",
        index = "star/{sample}.sorted.bam.bai"
    threads: 24
    # singularity:
    #     "shub://deanpettinga/rnaseq:rnaseq"
    log:
        "logs/sort_index_bam.{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools sort -@ {threads} -O BAM -o {output.sorted} {input}
        samtools index -b {output.sorted}
        '''
