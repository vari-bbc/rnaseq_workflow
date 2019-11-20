def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="trimmed_data/{sample}-SE_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "trimmed_data/{sample}-R1_val_1.fq.gz".format(**wildcards)
        fq2 = "trimmed_data/{sample}-R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule STAR:
    input:
        STAR_input
    output:
        # see STAR manual for additional output files
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
        "analysis/star/{sample}.Log.final.out",
        "analysis/star/{sample}.Log.out",
        "analysis/star/{sample}.Log.progress.out",
        "analysis/star/{sample}.ReadsPerGene.out.tab",
        "analysis/star/{sample}.SJ.out.tab",
        directory("analysis/star/{sample}._STARgenome"),
        directory("analysis/star/{sample}._STARpass1"),

    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts "
    threads: 24
    wrapper:
        "file:wrappers/star"
