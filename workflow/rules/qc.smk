def concat_fqs_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fqs = expand("results/trimmed_data/{sample}_{group_index}_R{read}_trimmed.fq.gz", group_index=units[units['sample'] == wildcards.sample]['group_index'], sample=wildcards.sample, read=wildcards.read)
    elif config["PE_or_SE"] == "PE":
        fqs =  expand("results/trimmed_data/{sample}_{group_index}_R{read}_val_{read}.fq.gz", group_index=units[units['sample'] == wildcards.sample]['group_index'], sample=wildcards.sample, read=wildcards.read)
    return fqs

rule concat_fastqs:
    """
    Concatenate fastqs.
    """
    input:
        concat_fqs_input
    output:
        "results/concat_fastqs/{sample}_R{read,[12]}.fastq.gz"
    benchmark:
        "benchmarks/concat_fastqs/{sample}_R{read}.txt"
    params:
        cat_or_symlink=lambda wildcards, input: "cat " + " ".join(input) + " > " if len(input) > 1 else "ln -sr " + input[0]
    threads: 1
    resources:
        mem_gb=4,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
    shell:
        """
        {params.cat_or_symlink} {output}
        """

rule fastq_screen:
    input:
                    "results/concat_fastqs/{fq_pref}.fastq.gz",
    output:
        html =      "results/fastq_screen/{fq_pref}_screen.html",
        txt =       "results/fastq_screen/{fq_pref}_screen.txt",
    benchmark:
                    "benchmarks/fastq_screen/{fq_pref}.bmk"
    threads: 8
    resources:
        nodes =     1,
        mem_gb =    32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules: config['modules']['fastq_screen']
    shell:
        """
        fastq_screen --threads {threads} --outdir results/fastq_screen/ {input} 
        """

rule fastqc:
    """
    Run fastqc on fastq files.
    """
    input:
        "results/concat_fastqs/{fq_pref}.fastq.gz"
    output:
        html="results/fastqc/{fq_pref}_fastqc.html",
        zip="results/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="results/fastqc/"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule seqtk:
    input:
        "results/concat_fastqs/{fq_pref}.fastq.gz"
    output:
        temp("results/subsample/{fq_pref}.fastq.gz"),
    envmodules:
        config['modules']['seqtk']
    params:
        num_subsamp = 50000,
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        seqtk sample -s 100 {input} {params.num_subsamp} | gzip -c > {output}
        """

def sortmerna_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="results/subsample/{sample}_R1.fastq.gz".format(**wildcards)
        return fq1
    if config["PE_or_SE"] == "PE":
        fq1 = "results/subsample/{sample}_R1.fastq.gz".format(**wildcards)
        fq2 = "results/subsample/{sample}_R2.fastq.gz".format(**wildcards)
        return [fq1,fq2]

rule sortmerna:
    input:
        sortmerna_input,
    output:
        directory("results/sortmerna/{sample}")
    envmodules:
        config['modules']['sortmerna']
    params:
        rfam5_8s = config["sortmerna"]["rfam5_8s"],
        rfam5s = config['sortmerna']['rfam5s'],
        silva_arc_16s = config['sortmerna']['silva_arc_16s'],
        silva_arc_23s = config['sortmerna']['silva_arc_23s'],
        silva_bac_16s = config['sortmerna']['silva_bac_16s'],
        silva_bac_23s = config['sortmerna']['silva_bac_23s'],
        silva_euk_18s = config['sortmerna']['silva_euk_18s'],
        silva_euk_28s = config['sortmerna']['silva_euk_28s'],
        idx_dir = config['sortmerna']['idx_dir'],
        fastqs = lambda wildcards, input: ' '.join(["-reads " + x for x in input])
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        sortmerna --threads {threads} {params.fastqs} --workdir {output}  \
        --idx-dir {params.idx_dir}  \
        --ref {params.rfam5s}  \
        --ref {params.rfam5_8s}  \
        --ref {params.silva_arc_16s}  \
        --ref {params.silva_arc_23s}  \
        --ref {params.silva_bac_16s}  \
        --ref {params.silva_bac_23s}  \
        --ref {params.silva_euk_18s}  \
        --ref {params.silva_euk_28s}
        """

multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R1_trimmed_screen.txt", samples=samples.itertuples()))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R1_trimmed_fastqc.zip", samples=samples.itertuples()))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R1_trimmed_fastqc.html", samples=samples.itertuples()))
    multiqc_input.append(expand("results/star/{samples.sample}.Log.final.out", samples=samples.itertuples()))
    multiqc_input.append(expand("results/sortmerna/{samples.sample}",samples=samples.itertuples()))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["aux_info/meta_info.json"]))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R{read}_screen.txt", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R{read}_fastqc.html", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R{read}_fastqc.zip", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/star/{samples.sample}.Log.final.out", samples=samples.itertuples()))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["libParams/flenDist.txt","aux_info/meta_info.json"]))
    multiqc_input.append(expand("results/sortmerna/{samples.sample}",samples=samples.itertuples()))


rule multiqc:
    input:
        multiqc_input
    params:
        "results/star/",
        "results/fastq_screen/",
        "results/fastqc/",
        "results/salmon/",
        "results/sortmerna/",
    output:
        "results/multiqc/multiqc_report.html",
        "results/multiqc/multiqc_report_data/multiqc.log",
        "results/multiqc/multiqc_report_data/multiqc_fastqc.txt",
        "results/multiqc/multiqc_report_data/multiqc_general_stats.txt",
        "results/multiqc/multiqc_report_data/multiqc_sources.txt",
        "results/multiqc/multiqc_report_data/multiqc_star.txt",
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 32,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    envmodules:
        config['modules']['multiqc']
    shell:
        """
        multiqc -f {params} \
        -o results/multiqc \
        --ignore '*._STARpass1/*' \
        -n multiqc_report.html \
        --cl-config 'max_table_rows: 999999' 
        """

