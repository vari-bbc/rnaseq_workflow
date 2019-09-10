rule sort_index_bam:
    input:
        bam="analysis/star/{sample}_{unit}.Aligned.out.bam",
        log="analysis/star/{sample}_{unit}.Log.final.out",
    output:
        sorted = "analysis/star/{sample}_{unit}.sorted.bam",
        index = "analysis/star/{sample}_{unit}.sorted.bam.bai"
    threads: 24
    conda:
        "../envs/samtools.yaml"
    shell:
        '''
        samtools sort -@ {threads} -O BAM -o {output.sorted} {input.bam}
        samtools index -b {output.sorted}
        '''
