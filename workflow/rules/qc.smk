rule fastq_screen:
    input:
                    "results/rename_fastqs/{fq_pref}.fastq.gz",
    output:
        html =      "results/fastq_screen/{fq_pref}_screen.html",
        txt =       "results/fastq_screen/{fq_pref}_screen.txt",
    benchmark:
                    "benchmarks/fastq_screen/{fq_pref}.bmk"
    threads: 8
    resources:
        nodes =     1,
        mem_gb =    32,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules: config['modules']['fastq_screen']
    shell:
        """
        fastq_screen --threads {threads} --outdir results/fastq_screen/ {input} 
        """

rule seqtk_SE:
    input:
        fq1 = STAR_input,
    output:
        fq1 = temp("results/subsample/{sample}_R1_trimmed.fq.gz"),
    envmodules:
        config['modules']['seqtk']
    params:
        num_subsamp = 50000,
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        seqtk sample -s 100 {input.fq1} {params.num_subsamp} | gzip -c > {output.fq1}
        """

rule seqtk_PE:
    input:
        fq1 = lambda wildcards:  STAR_input(wildcards)[0],
        fq2 = lambda wildcards:  STAR_input(wildcards)[1],
    output:
        fq1 = temp("results/subsample/{sample}_R1_val_1.fq.gz"),
        fq2 = temp("results/subsample/{sample}_R2_val_2.fq.gz"),
    envmodules:
        config['modules']['seqtk']
    params:
        num_subsamp = 50000,
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        seqtk sample -s 100 {input.fq1} {params.num_subsamp} | gzip -c > {output.fq1}
        seqtk sample -s 100 {input.fq2} {params.num_subsamp} | gzip -c > {output.fq2}
        """

def sortmerna_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="results/subsample/{sample}_R1_trimmed.fq.gz".format(**wildcards)
        return fq1
    if config["PE_or_SE"] == "PE":
        fq1 = "results/subsample/{sample}_R1_val_1.fq.gz".format(**wildcards)
        fq2 = "results/subsample/{sample}_R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule sortmerna_SE:
    input:
        fq1 = sortmerna_input,
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
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        sortmerna --threads {threads} -reads {input.fq1} --workdir {output}  \
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

rule sortmerna_PE:
    input:
        fq1 = lambda wildcards: sortmerna_input(wildcards)[0],
        fq2 = lambda wildcards: sortmerna_input(wildcards)[1],
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
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        sortmerna --threads {threads} -reads {input.fq1} -reads {input.fq2} --workdir {output}  \
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
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R1_screen.txt", samples=samples.itertuples()))
    multiqc_input.append(expand("results/trimmed_data/{samples.sample}_R1_trimmed_fastqc.zip", samples=samples.itertuples()))
    multiqc_input.append(expand("results/trimmed_data/{samples.sample}_R1_trimmed_fastqc.html", samples=samples.itertuples()))
    multiqc_input.append(expand("results/star/{samples.sample}.Log.final.out", samples=samples.itertuples()))
    multiqc_input.append(expand("results/sortmerna/{samples.sample}",samples=samples.itertuples()))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["aux_info/meta_info.json"]))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R{read}_screen.txt", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/trimmed_data/{samples.sample}_R{read}_val_{read}_fastqc.html", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/trimmed_data/{samples.sample}_R{read}_val_{read}_fastqc.zip", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/trimmed_data/{samples.sample}_R{read}.fastq.gz_trimming_report.txt", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/star/{samples.sample}.Log.final.out", samples=samples.itertuples()))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["libParams/flenDist.txt","aux_info/meta_info.json"]))
    multiqc_input.append(expand("results/sortmerna/{samples.sample}",samples=samples.itertuples()))


rule multiqc:
    input:
        multiqc_input
    params:
        "results/star/",
        "results/trimmed_data/",
        "results/fastq_screen/",
        "results/salmon/",
        "results/sortmerna/",
    output:
        "results/multiqc/multiqc_report.html",
        "results/multiqc/multiqc_report_data/multiqc.log",
        "results/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
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
        log_prefix=lambda wildcards: "_".join(wildcards)
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

