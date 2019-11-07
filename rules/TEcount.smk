rule TEcount:
    input:
        BAM="analysis/star/{sample}.Aligned.out.bam",
        GTF=config["ref"]["annotation"],
        TE_GTF=config["TEref"]["annotation"],

    output:
        "analysis/TEcount/{sample}.cntTable",
    params:
        strand=config["strandedness"],
        project="analysis/TEcount/{sample}",
        multi_or_uniq="uniq",
    conda:
        "../envs/TEtranscripts.yaml"
    shell:
        """
        TEcount \
        --BAM {input.BAM} \
        --GTF {input.GTF} \
        --TE {input.TE_GTF} \
        --format BAM \
        --stranded {params.strand} \
        --project {params.project} \
        --mode {params.multi_or_uniq}
        """
