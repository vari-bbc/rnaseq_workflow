def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq2"].values)
    return fastq

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "results/rename_fastqs/{sample}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
        cat_or_symlink=lambda wildcards, input: "cat " + " ".join(input) + " > " if len(input) > 1 else "ln -sr " + input[0]
    threads: 1
    resources:
        mem_gb=4,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        {params.cat_or_symlink} {output}
        """

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "results/rename_fastqs/{sample}_R1.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 = "results/rename_fastqs/{sample}_R1.fastq.gz".format(**wildcards)
        R2 = "results/rename_fastqs/{sample}_R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        temp("results/trimmed_data/{sample}_R1_val_1.fq.gz"),
        "results/trimmed_data/{sample}_R1_val_1_fastqc.html",
        "results/trimmed_data/{sample}_R1_val_1_fastqc.zip",
        "results/trimmed_data/{sample}_R1.fastq.gz_trimming_report.txt",
        temp("results/trimmed_data/{sample}_R2_val_2.fq.gz"),
        "results/trimmed_data/{sample}_R2_val_2_fastqc.html",
        "results/trimmed_data/{sample}_R2_val_2_fastqc.zip",
        "results/trimmed_data/{sample}_R2.fastq.gz_trimming_report.txt"
    params:
        outdir="results/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        trim_galore \
        --paired \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q 20 \
        --fastqc
        """

rule trim_galore_SE:
    input:
        trim_galore_input
    output:
        temp("results/trimmed_data/{sample}_R1_trimmed.fq.gz"),
        "results/trimmed_data/{sample}_R1_trimmed_fastqc.zip",
        "results/trimmed_data/{sample}_R1_trimmed_fastqc.html",
        "results/trimmed_data/{sample}_R1.fastq.gz_trimming_report.txt",
    params:
        outdir="results/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        trim_galore \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q 20 \
        --fastqc 
        """

def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="results/trimmed_data/{sample}_R1_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "results/trimmed_data/{sample}_R1_val_1.fq.gz".format(**wildcards)
        fq2 = "results/trimmed_data/{sample}_R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule STAR:
    input:
        STAR_input
    output:
        # see STAR manual for additional output files
        bam =                 "results/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai =                 "results/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final =           "results/star/{sample}.Log.final.out",
        log =                 "results/star/{sample}.Log.out",
        rpg =                 "results/star/{sample}.ReadsPerGene.out.tab",
        sj =                  "results/star/{sample}.SJ.out.tab",
        g_dir =     directory("results/star/{sample}._STARgenome"),
        pass1_dir = directory("results/star/{sample}._STARpass1"),
    params:
        # path to STAR reference genome index
        index = config["ref"]["index"],
        outprefix = "results/star/{sample}."
    benchmark:
        "benchmarks/star/{sample}.txt"
    envmodules:
        config['modules']['star'],
        config['modules']['samtools']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {input} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log 

        samtools index {output.bam}
        """

rule salmon:
    input:
        STAR_input,
        index=config["ref"]["salmon_index"]
    output:
        expand("results/salmon/{{sample}}/{file}", file=["libParams/flenDist.txt","aux_info/meta_info.json","quant.sf","lib_format_counts.json","cmd_info.json","logs/salmon_quant.log"])
    params:
        outdir=directory("results/salmon/{sample}"),
        reads=lambda wildcards, input: "-1 {fq1} -2 {fq2}".format(fq1=input[0], fq2=input[1]) if config["PE_or_SE"] == "PE" else "-r {fq1}".format(fq1=input[0])
    benchmark:
        "benchmarks/salmon/{sample}.txt"
    envmodules:
        config['modules']['salmon']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        salmon quant \
                -p {threads} \
                -l A \
                -i {input.index} \
                {params.reads} \
                --validateMappings \
                -o {params.outdir}
        """

