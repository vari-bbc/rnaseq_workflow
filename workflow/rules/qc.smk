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

rule make_genes_ref_flat:
    """
    Make REF FLAT file from GTF.
    """
    input:
        gtf=config['ref']['annotation']
    output:
        ref_flat="results/misc/gene_annot.ref_flat.txt"
    benchmark:
        "benchmarks/make_genes_ref_flat/bench.txt"
    params:
    envmodules:
        config["modules"]["ucsctools"]
    threads: 4
    resources:
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons {input.gtf} /dev/stdout | \
                perl -F"\\t" -lane 'print join("\\t", @F[11,0..9])'  > {output.ref_flat}
        """

rule make_genes_bed:
    """
    Make BED file from GTF.
    """
    input:
        gtf=config['ref']['annotation']
    output:
        bed="results/misc/gene_annot.bed"
    benchmark:
        "benchmarks/make_genes_bed/bench.txt"
    params:
    envmodules:
        config["modules"]["ucsctools"]
    threads: 4
    resources:
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gtfToGenePred {input.gtf} /dev/stdout |  genePredToBed /dev/stdin {output.bed} 
        """


rule get_rRNA_intervals_from_gtf:
    """
    Get rRNA intervals from GTF in BED format.
    """
    input:
        gtf=config['ref']['annotation'],
        ref_dict=config['ref']['dict'],
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname'])),
    output:
        bed="results/misc/rrna.bed",
        interval_list="results/misc/rrna.interval_list"
    benchmark:
        "benchmarks/get_rRNA_intervals_from_gtf/bench.txt"
    params:
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock),
    envmodules:
        config["modules"]["picard"],
        config["modules"]["R"],
    threads: 4
    resources:
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        Rscript --vanilla -e 'renv::load("{params.renv_rproj_dir}"); library(rtracklayer); gtf <- import("{input.gtf}"); biotype_col <- na.omit(match(c("gene_biotype","gene_type"), colnames(mcols(gtf)))); stopifnot(length(biotype_col) > 0); rrna <- gtf[mcols(gtf)[[biotype_col[1]]]=="rRNA" & mcols(gtf)$type=="gene"]; score(rrna) <- 1; export(rrna, "{output.bed}")'

        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD BedToIntervalList \
                I={output.bed} \
                O={output.interval_list} \
                SD={input.ref_dict}
        """

def get_library_strandedness(wildcards, input):
    with open(input.strandedness) as f:
        first_line = f.readline().strip('\n')
    if first_line in ("ISR","SR"):
        strandedness = "SECOND_READ_TRANSCRIPTION_STRAND"
    elif first_line in ("ISF","SF"):
        strandedness = "FIRST_READ_TRANSCRIPTION_STRAND"
    elif first_line in ("IU","U"):
        strandedness = "NONE"
    else:
        raise Exception("Unrecognized strand type") 
    return strandedness

rule CollectRnaSeqMetrics:
    """
    Run Picard CollectRnaSeqMetrics.
    """
    input:
        bam="results/star/{sample}.sorted.bam",
        ref_flat="results/misc/gene_annot.ref_flat.txt",
        strandedness="results/SummarizedExperiment/inferred_strandedness.txt",
        rrna="results/misc/rrna.interval_list"
    output:
        metrics="results/CollectRnaSeqMetrics/{sample}.txt"
    benchmark:
        "benchmarks/CollectRnaSeqMetrics/{sample}.txt"
    params:
        strand=get_library_strandedness,
    envmodules:
        config["modules"]["picard"]
    threads: 4
    resources:
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD CollectRnaSeqMetrics \
        -I {input.bam} \
        -O {output.metrics} \
        --REF_FLAT {input.ref_flat} \
        --STRAND_SPECIFICITY {params.strand} \
        --RIBOSOMAL_INTERVALS {input.rrna}
        """

rule rseqc_genebody_cov:
    """
    Run RSeQC.
    """
    input:
        bam="results/star/{sample}.sorted.bam",
        genes="results/misc/gene_annot.bed"
    output:
        metrics="results/rseqc_genebody_cov/{sample}/{sample}.geneBodyCoverage.txt",
        r="results/rseqc_genebody_cov/{sample}/{sample}.geneBodyCoverage.r",
        log="results/rseqc_genebody_cov/{sample}/log.txt",
    benchmark:
        "benchmarks/rseqc_genebody_cov/{sample}.txt"
    params:
        prefix="{sample}",
        samp_dir=lambda wildcards, output: os.path.dirname(output.metrics)
    envmodules:
        config["modules"]["rseqc"]
    threads: 4
    resources:
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        cd {params.samp_dir}

        geneBody_coverage.py -r ../../../{input.genes} -i ../../../{input.bam} -o {params.prefix}
        """

multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R1_trimmed_screen.txt", samples=samples.itertuples()))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R1_trimmed_fastqc.zip", samples=samples.itertuples()))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R1_trimmed_fastqc.html", samples=samples.itertuples()))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["aux_info/meta_info.json"]))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("results/fastq_screen/{samples.sample}_R{read}_screen.txt", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R{read}_fastqc.html", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/fastqc/{samples.sample}_R{read}_fastqc.zip", samples=samples.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("results/salmon/{samples.sample}/{file}", samples=samples.itertuples(), file=["libParams/flenDist.txt","aux_info/meta_info.json"]))

multiqc_input.append(expand("results/star/{samples.sample}.Log.final.out", samples=samples.itertuples()))
multiqc_input.append(expand("results/sortmerna/{samples.sample}",samples=samples.itertuples()))
multiqc_input.append(expand("results/CollectRnaSeqMetrics/{samples.sample}.txt",samples=samples.itertuples()))
if config['run_rseqc']:
    multiqc_input.append(expand("results/rseqc_genebody_cov/{samples.sample}/{samples.sample}.geneBodyCoverage.txt",samples=samples.itertuples()))


rule multiqc:
    input:
        multiqc_input
    params:
        lambda wildcards, input: " ".join(pd.unique([os.path.dirname(x) for x in input]))
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

