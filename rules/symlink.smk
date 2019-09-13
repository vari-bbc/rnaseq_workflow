def symlink_input(wildcards):
    if config["PE_or_SE"]=="SE":
        SE_reads = units.loc[:,"fq1"].values.tolist()
        SE_reads = ["raw_reads/"+read for read in SE_reads]
        return SE_reads
    elif config["PE_or_SE"]=="PE":
        PE_reads = units.loc[:,["fq1","fq2"]].values.tolist()
        # flatten the list of lists
        PE_reads = [read for pair in PE_reads for read in pair]
        PE_reads = ["raw_reads/"+read for read in PE_reads]
        return PE_reads

def get_fq(wildcards):
    if config["PE_or_SE"]=="SE":
        reads = units.loc[lambda df: df['sample'] == wildcards.sample].loc[lambda df: df['unit'] == wildcards.unit,"fq1"].values.tolist()[0]
        reads = "raw_reads/"+reads
    elif config["PE_or_SE"]=="PE":
        reads = units.loc[lambda df: df['sample'] == wildcards.sample].loc[lambda df: df['unit'] == wildcards.unit,("fq1","fq2")].values.tolist()[0]
        for sample in reads:
            sample = "raw_reads/"+sample
    return reads

rule symlink_SE:
    # input:
    #     symlink_input
    output:
        expand("raw_reads/{sample}-{unit}.fastq.gz",sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
        # expand("raw_reads/{sample}_{unit}_R1.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
        # expand("raw_reads/{sample}_{unit}_R2.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
    priority: 2
    shell:
        "Rscript src/symlink_SE.R"

rule symlink_PE:
    # input:
    #     symlink_input
    output:
        expand("raw_reads/{sample}_{unit}_R1.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
        expand("raw_reads/{sample}_{unit}_R2.fastq.gz", sample=set(units["sample"].tolist()), unit=set(units["unit"].tolist())),
    priority: 2
    shell:
        "Rscript src/symlink.R"
