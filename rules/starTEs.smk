def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="trimmed_data/{sample}-SE_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "trimmed_data/{sample}-R1_val_1.fq.gz".format(**wildcards)
        fq2 = "trimmed_data/{sample}-R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule starTEs:
    input:
        STAR_input,
        "ref/gencode/Genome",
        "ref/gencode/Log.out",
        "ref/gencode/SA",
        "ref/gencode/SAindex",
        "ref/gencode/chrLength.txt",
        "ref/gencode/chrName.txt",
        "ref/gencode/chrNameLength.txt",
        "ref/gencode/chrStart.txt",
        "ref/gencode/exonGeTrInfo.tab",
        "ref/gencode/exonInfo.tab",
        "ref/gencode/geneInfo.tab",
        "ref/gencode/genomeParameters.txt",
        "ref/gencode/sjdbInfo.txt",
        "ref/gencode/sjdbList.fromGTF.out.tab",
        "ref/gencode/sjdbList.out.tab",
        "ref/gencode/transcriptInfo.tab",
    output:
        # see STAR manual for additional output files
        "analysis/starTE/{sample}.Aligned.out.bam",
        "analysis/starTE/{sample}.Log.final.out",
        "analysis/starTE/{sample}.Log.out",
        "analysis/starTE/{sample}.Log.progress.out",
        "analysis/starTE/{sample}.ReadsPerGene.out.tab",
        "analysis/starTE/{sample}.SJ.out.tab",
        directory("analysis/starTE/{sample}._STARgenome"),
        # directory("analysis/star/{sample}._STARpass1"),

    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        GTF=config["ref"]["annotation"],
        readLength=config["readLength"],
        # optional parameters
        # extra="--quantMode GeneCounts "
    threads: 24
    wrapper:
        "file:wrappers/starTEs/wrapper.py"
