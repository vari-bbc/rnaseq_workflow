rule align:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = "trimmed_data/{sample}_{unit}_R1_val_1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "trimmed_data/{sample}_{unit}_R2_val_2.fq.gz"
    output:
        # see STAR manual for additional output files
        "analysis/star/{sample}_{unit}.ReadsPerGene.out.tab",
        "analysis/star/{sample}_{unit}.Aligned.out.bam",
        "analysis/star/{sample}_{unit}.Log.out",
        "analysis/star/{sample}_{unit}.Log.progress.out",
        "analysis/star/{sample}_{unit}.Log.final.out",
    log:
        "logs/star/{sample}_{unit}.log"
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
        "analysis/star/{sample}_{unit}.Aligned.out.bam",
    output:
        sorted = "analysis/star/{sample}_{unit}.sorted.bam",
        index = "analysis/star/{sample}_{unit}.sorted.bam.bai"
    threads: 24
    # singularity:
    #     "shub://deanpettinga/rnaseq:rnaseq"
    log:
        "logs/sort/sort_index_bam.{sample}_{unit}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools sort -@ {threads} -O BAM -o {output.sorted} {input}
        samtools index -b {output.sorted}
        '''
