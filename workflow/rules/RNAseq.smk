import gzip

def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = units[(units["sample"] == wildcards.sample) & (units["group_index"] == wildcards.group_index)]["fq1"]
    elif wildcards.read == "R2":
            fastq = units[(units["sample"] == wildcards.sample) & (units["group_index"] == wildcards.group_index)]["fq2"]
    return 'raw_data/' + fastq

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name.
    """
    input:
        get_orig_fastq
    output:
        "results/rename_fastqs/{sample}_{group_index}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{group_index}_{read}.txt"
    params:
    threads: 1
    resources:
        mem_gb=4,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
    shell:
        """
        ln -sr {input} {output} 
        """

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "results/rename_fastqs/{sample}_{group_index}_R1.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 = "results/rename_fastqs/{sample}_{group_index}_R1.fastq.gz".format(**wildcards)
        R2 = "results/rename_fastqs/{sample}_{group_index}_R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        #temp("results/trimmed_data/{sample}_{group_index}_R1_val_1.fq.gz"),
        "results/trimmed_data/{sample}_{group_index}_R1_val_1.fq.gz",
        "results/trimmed_data/{sample}_{group_index}_R1_val_1_fastqc.html",
        "results/trimmed_data/{sample}_{group_index}_R1_val_1_fastqc.zip",
        "results/trimmed_data/{sample}_{group_index}_R1.fastq.gz_trimming_report.txt",
        #temp("results/trimmed_data/{sample}_{group_index}_R2_val_2.fq.gz"),
        "results/trimmed_data/{sample}_{group_index}_R2_val_2.fq.gz",
        "results/trimmed_data/{sample}_{group_index}_R2_val_2_fastqc.html",
        "results/trimmed_data/{sample}_{group_index}_R2_val_2_fastqc.zip",
        "results/trimmed_data/{sample}_{group_index}_R2.fastq.gz_trimming_report.txt"
    params:
        outdir="results/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}_{group_index}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
        temp("results/trimmed_data/{sample}_{group_index}_R1_trimmed.fq.gz"),
        "results/trimmed_data/{sample}_{group_index}_R1_trimmed_fastqc.zip",
        "results/trimmed_data/{sample}_{group_index}_R1_trimmed_fastqc.html",
        "results/trimmed_data/{sample}_{group_index}_R1.fastq.gz_trimming_report.txt",
    params:
        outdir="results/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}_{group_index}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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
    group_indices = units[units['sample'] == wildcards.sample]['group_index']
    if config["PE_or_SE"] == "SE":
        fq1 = expand("results/trimmed_data/{sample}_{group_index}_R1_trimmed.fq.gz", sample=wildcards.sample, group_index=group_indices)
        fq2 = []
    elif config["PE_or_SE"] == "PE":
        fq1 = expand("results/trimmed_data/{sample}_{group_index}_R1_val_1.fq.gz", sample=wildcards.sample, group_index=group_indices)
        fq2 = expand("results/trimmed_data/{sample}_{group_index}_R2_val_2.fq.gz", sample=wildcards.sample, group_index=group_indices)
    return {'fq1_files': fq1, 'fq2_files': fq2}

def get_RG(wildcards, input):
    fq1_files = input.fq1_files

    ## Use the user-specified read group info if available
    rg_lines = units[units['sample'] == wildcards.sample]['RG'].values
    
    if(pd.isnull(rg_lines).any()):

        rg_lines = []
        
        ## Extract the first line of each fq1 file
        first_lines = []
        for fq_file in fq1_files:
            with gzip.open(fq_file,'rt') as f:
                first_lines.append(f.readline().strip())

        ## Compile the read group line for each library
        for i in range(len(fq1_files)):
            first_line_split = first_lines[i].split(':')
         
            flowcell = first_line_split[2]
            lane = first_line_split[3]
            lib_barcode = first_line_split[9]
         
            sample = wildcards.sample
         
            rgid = '.'.join([flowcell, lane, lib_barcode])
            rgpu = '.'.join([flowcell, lane, lib_barcode])
            rgsm = sample

            # Assumes one library per sample
            rg_line = 'ID:' + rgid + ' PU:' + rgpu + ' LB:' + rgsm + ' PL:ILLUMINA SM:' + rgsm
            rg_lines.append(rg_line)

    read_groups = ' , '.join(rg_lines)

    return read_groups

rule STAR:
    input:
        unpack(STAR_input)
    output:
        # see STAR manual for additional output files
        unsorted_bam =        temp("results/star/{sample}.Aligned.out.bam"),
        log_final =           "results/star/{sample}.Log.final.out",
        log =                 "results/star/{sample}.Log.out",
        rpg =                 "results/star/{sample}.ReadsPerGene.out.tab",
        sj =                  "results/star/{sample}.SJ.out.tab",
        g_dir =     directory("results/star/{sample}._STARgenome"),
        pass1_dir = directory("results/star/{sample}._STARpass1"),
        bam =                 "results/star/{sample}.sorted.bam",
        bai =                 "results/star/{sample}.sorted.bam.bai",
    params:
        # path to STAR reference genome index
        index = config["ref"]["index"],
        outprefix = "results/star/{sample}.",
        in_fastqs = lambda wildcards, input: ','.join(input.fq1_files) + ' ' + ','.join(input.fq2_files),
        read_groups = get_RG
    benchmark:
        "benchmarks/star/{sample}.txt"
    envmodules:
        config['modules']['star'],
        config['modules']['samtools']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.in_fastqs} \
        --outSAMattrRGline {params.read_groups} \
        --twopassMode Basic \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix} \
        --quantMode GeneCounts \
        --outStd Log 

        samtools sort -@ {threads} -o {output.bam} {output.unsorted_bam}
        samtools index {output.bam}
        """

rule salmon:
    input:
        unpack(STAR_input),
        index=config["ref"]["salmon_index"]
    output:
        expand("results/salmon/{{sample}}/{file}", file=["libParams/flenDist.txt","aux_info/meta_info.json","quant.sf","lib_format_counts.json","cmd_info.json","logs/salmon_quant.log"])
    params:
        outdir=directory("results/salmon/{sample}"),
        reads=lambda wildcards, input: "-1 {fq1} -2 {fq2}".format(fq1=input.fq1_files, fq2=input.fq2_files) if config["PE_or_SE"] == "PE" else "-r {fq1}".format(fq1=input.fq1_files)
    benchmark:
        "benchmarks/salmon/{sample}.txt"
    envmodules:
        config['modules']['salmon']
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
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

rule SummarizedExperiment:
    input:
        star = expand("results/star/{samples.sample}.ReadsPerGene.out.tab", samples=samples.itertuples()),
        salmon = expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=['lib_format_counts.json', 'quant.sf'])
    output:
        se="results/SummarizedExperiment/SummarizedExperiment.rds",
        sce="results/SummarizedExperiment/sce.rds",
        sizeFactors="results/SummarizedExperiment/DESeq2_sizeFactors_reciprocal.tsv"
    benchmark:
        "benchmarks/SummarizedExperiment/SummarizedExperiment.txt"
    params:
        gtf=config['ref']['annotation'],
        orgdb=config['orgdb']
    threads: 1
    envmodules:
        config['modules']['R']
    resources:
        mem_gb=64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        Rscript --vanilla workflow/scripts/make_sce.R {params.gtf} {params.orgdb} {output.se} {output.sce} {output.sizeFactors}
        """
