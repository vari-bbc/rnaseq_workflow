def get_fq(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="trimmed_data/{sample}-{unit}_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "trimmed_data/{sample}_{unit}_R1_val_1.fq.gz".format(**wildcards)
        fq2 = "trimmed_data/{sample}_{unit}_R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule star_align:
    input:
        get_fq
    output:
        # see STAR manual for additional output files
        "analysis/star/{sample}_{unit}.Aligned.out.bam",
        "analysis/star/{sample}_{unit}.Log.final.out",
        "analysis/star/{sample}_{unit}.Log.out",
        "analysis/star/{sample}_{unit}.Log.progress.out",
        "analysis/star/{sample}_{unit}.ReadsPerGene.out.tab",
        "analysis/star/{sample}_{unit}.SJ.out.tab",
        directory("analysis/star/{sample}_{unit}._STARgenome"),
        directory("analysis/star/{sample}_{unit}._STARpass1"),

    log:
        "logs/star/{sample}_{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts "
    threads: 24
    wrapper:
        "file:wrappers/star/wrapper.py"
