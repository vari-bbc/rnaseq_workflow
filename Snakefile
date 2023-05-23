import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("7.25.0")

##### check if conda is available. If not, we don't run rules requiring conda envs. Applies to dryruns and actual runs.
conda_avail = (which('conda') is not None)

##### load config and sample sheets #####

configfile: "bin/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

# units = pd.read_table("bin/units_test.tsv").set_index("sample", drop=False)
units = pd.read_table(config["units"]).set_index("sample", drop=False)
var_calling_units = pd.read_table("bin/variant_calling_units.tsv").set_index("unit", drop=False)


if not os.path.exists('tmp'):
    os.mkdir('tmp')

fai_file = config["ref"]["fai"]
contigs_file = "bin/grouped_contigs.tsv"

# read in file with col1 as contig group name and col2 as a commma separated list of contigs. This was written for use with GATK for splitting the variant calling steps by chromosome/contig. We group the unplaced contigs together since those tend to be small.
# we read in this file even if not calling variants. Otherwise variant-calling rules relying on this file will cause an error.
contig_groups = pd.read_table(contigs_file)

# if config set to run variants component of pipeline,  check that chromosomes were parsed correctly in the grouped_Contigs.tsv file. 
if config["call_variants"]:

    # check chromosomes/contigs parsed correctly.
    chroms = contig_groups['contigs'].values[0:(contig_groups.shape[0]-1)].tolist()
    unanchored_contigs = contig_groups['contigs'].values[contig_groups.shape[0]-1].split(",")
    contig_groups_flattened = chroms + unanchored_contigs
    
    fai_pd = pd.read_table(fai_file, header=None)

    assert contig_groups_flattened == fai_pd.iloc[:,0].values.tolist(), "Chromosomes in .fai do not match those in 'grouped_contigs.tsv'."

avail_read = ["1", "2"] if config["PE_or_SE"] == "PE" else ["1"]

##### target rules #####

# Need this directive because both PE and SE versions of these rules produce the trimmed R1 output files.
ruleorder: trim_galore_PE > trim_galore_SE
ruleorder: seqtk_PE > seqtk_SE
ruleorder: sortmerna_PE > sortmerna_SE

rule all:
    input:
        lambda wildcards: "analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi" if config["call_variants"] else [],
        "analysis/09b_snp_pca_and_dendro/report.html" if (config["call_variants"] and conda_avail) else [],
        "analysis/multiqc/multiqc_report.html",


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
        "analysis/rename_fastqs/{sample}_{read}.fastq.gz"
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
        cat_or_symlink=lambda wildcards, input: "cat " + " ".join(input) + " > " if len(input) > 1 else "ln -sr " + input[0]
    threads: 1
    resources:
        mem_gb=8,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
    shell:
        """
        {params.cat_or_symlink} {output}
        """

rule fastq_screen:
    input:
                    "analysis/rename_fastqs/{fq_pref}.fastq.gz",
    output:
        html =      "analysis/fastq_screen/{fq_pref}_screen.html",
        txt =       "analysis/fastq_screen/{fq_pref}_screen.txt",
    benchmark:
                    "benchmarks/fastq_screen/{fq_pref}.bmk"
    threads: 8
    resources:
        nodes =     1,
        mem_gb =    64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules: config['modules']['fastq_screen']
    shell:
        """
        fastq_screen --outdir analysis/fastq_screen/ {input} 
        """

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "analysis/rename_fastqs/{sample}_R1.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 = "analysis/rename_fastqs/{sample}_R1.fastq.gz".format(**wildcards)
        R2 = "analysis/rename_fastqs/{sample}_R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        "analysis/trimmed_data/{sample}_R1_val_1.fq.gz",
        "analysis/trimmed_data/{sample}_R1_val_1_fastqc.html",
        "analysis/trimmed_data/{sample}_R1_val_1_fastqc.zip",
        "analysis/trimmed_data/{sample}_R1.fastq.gz_trimming_report.txt",
        "analysis/trimmed_data/{sample}_R2_val_2.fq.gz",
        "analysis/trimmed_data/{sample}_R2_val_2_fastqc.html",
        "analysis/trimmed_data/{sample}_R2_val_2_fastqc.zip",
        "analysis/trimmed_data/{sample}_R2.fastq.gz_trimming_report.txt"
    params:
        outdir="analysis/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  80,
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
        "analysis/trimmed_data/{sample}_R1_trimmed.fq.gz",
        "analysis/trimmed_data/{sample}_R1_trimmed_fastqc.zip",
        "analysis/trimmed_data/{sample}_R1_trimmed_fastqc.html",
        "analysis/trimmed_data/{sample}_R1.fastq.gz_trimming_report.txt",
    params:
        outdir="analysis/trimmed_data/"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  80,
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
        fq1="analysis/trimmed_data/{sample}_R1_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "analysis/trimmed_data/{sample}_R1_val_1.fq.gz".format(**wildcards)
        fq2 = "analysis/trimmed_data/{sample}_R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule STAR:
    input:
        STAR_input
    output:
        # see STAR manual for additional output files
        bam =                 "analysis/star/{sample}.Aligned.sortedByCoord.out.bam",
        bai =                 "analysis/star/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log_final =           "analysis/star/{sample}.Log.final.out",
        log =                 "analysis/star/{sample}.Log.out",
        rpg =                 "analysis/star/{sample}.ReadsPerGene.out.tab",
        sj =                  "analysis/star/{sample}.SJ.out.tab",
        g_dir =     directory("analysis/star/{sample}._STARgenome"),
        pass1_dir = directory("analysis/star/{sample}._STARpass1"),
    params:
        # path to STAR reference genome index
        index = config["ref"]["index"],
        outprefix = "analysis/star/{sample}."
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
        expand("analysis/salmon/{{sample}}/{file}", file=["libParams/flenDist.txt","aux_info/meta_info.json","quant.sf","lib_format_counts.json","cmd_info.json","logs/salmon_quant.log"])
    params:
        outdir=directory("analysis/salmon/{sample}"),
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


multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("analysis/fastq_screen/{units.sample}_R1_screen.txt", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R1_trimmed.fq.gz", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R1_trimmed_fastqc.zip", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R1_trimmed_fastqc.html", units=units.itertuples()))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
    multiqc_input.append(expand("analysis/sortmerna/{units.sample}",units=units.itertuples()))
    multiqc_input.append(expand("analysis/salmon/{units.sample}/{file}", units=units.itertuples(), file=["aux_info/meta_info.json"]))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("analysis/fastq_screen/{units.sample}_R{read}_screen.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R{read}_val_{read}.fq.gz", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R{read}_val_{read}_fastqc.html", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R{read}_val_{read}_fastqc.zip", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}_R{read}.fastq.gz_trimming_report.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
    multiqc_input.append(expand("analysis/salmon/{units.sample}/{file}", units=units.itertuples(), file=["libParams/flenDist.txt","aux_info/meta_info.json"]))
    multiqc_input.append(expand("analysis/sortmerna/{units.sample}",units=units.itertuples()))

rule multiqc:
    input:
        multiqc_input
    params:
        "analysis/star/",
        "analysis/trimmed_data/",
        "analysis/fastq_screen/",
        "analysis/salmon/",
        "analysis/sortmerna/",
    output:
        "analysis/multiqc/multiqc_report.html",
        "analysis/multiqc/multiqc_report_data/multiqc.log",
        "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_star.txt",
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
        -o analysis/multiqc \
        --ignore '*._STARpass1/*' \
        -n multiqc_report.html \
        --cl-config 'max_table_rows: 999999' 
        """

rule seqtk_SE:
    input:
        fq1 = STAR_input,
    output:
        fq1 = "analysis/subsample/{sample}_R1_trimmed.fq.gz",
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
        fq1 = "analysis/subsample/{sample}_R1_val_1.fq.gz",
        fq2 = "analysis/subsample/{sample}_R2_val_2.fq.gz",
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
        fq1="analysis/subsample/{sample}_R1_trimmed.fq.gz".format(**wildcards)
        return fq1
    if config["PE_or_SE"] == "PE":
        fq1 = "analysis/subsample/{sample}_R1_val_1.fq.gz".format(**wildcards)
        fq2 = "analysis/subsample/{sample}_R2_val_2.fq.gz".format(**wildcards)
        return [fq1,fq2]

rule sortmerna_SE:
    input:
        fq1 = sortmerna_input,
    output:
        directory("analysis/sortmerna/{sample}")
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
        directory("analysis/sortmerna/{sample}")
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


rule addRG:
    input:
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam="analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam",
    benchmark:
        "benchmarks/00_addRG/{sample}.txt"
    envmodules:
        config['modules']['picard']
    params:
        unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['unit'].values[0],
        lib=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['library'].values[0],
        sample=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['sample'].values[0],
        platform_unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['platform_unit'].values[0]
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms16g -Xmx64g -Djava.io.tmpdir=./tmp -jar $PICARD AddOrReplaceReadGroups \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --RGID {params.unit} \
        --RGLB {params.lib} \
        --RGPL ILLUMINA \
        --RGPU {params.platform_unit} \
        --RGSM {params.sample} 
        """


# variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020

rule markdups:
    input:
        lambda wildcards: expand("analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
    output:
        bam="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam",
        metrics="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.metrics"
    params: 
        input_params=lambda wildcards: expand("--INPUT analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
    benchmark:
        "benchmarks/01_markdup/{sample}.txt"
    envmodules:
        config['modules']['picard']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms16g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates \
        {params.input_params} \
        --OUTPUT {output.bam} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --METRICS_FILE {output.metrics} 
        """

rule splitncigar:
    input:
        "analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam"
    output:
        "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
    benchmark:
        "benchmarks/02_splitncigar/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input} \
        -O {output} 

        """

rule base_recalibrate:
    input:
        "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
    output:
        "analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
    benchmark:
        "benchmarks/03_base_recal/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        known_snps=config["ref"]["known_snps"],
        known_indels=config["ref"]["known_indels"]
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            BaseRecalibrator \
            -R {params.ref_fasta} \
            -I {input} \
            -O {output} \
            -known-sites {params.known_snps} \
            -known-sites {params.known_indels} 

        """

rule applyBQSR:
    input:
        bam="analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam",
        recal_table="analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
    output:
        "analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
    benchmark:
        "benchmarks/04_apply_base_recal/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 64,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {params.ref_fasta} \
            -I {input.bam} \
            -O {output} \
            --bqsr-recal-file {input.recal_table} 

        """

rule haplotypecaller:
    input:
        bam="analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
    output:
        "analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz"
    benchmark:
        "benchmarks/05_haplotypecaller/{sample}.{contig_group}.txt"
    params:
        dbsnp=config["ref"]["known_snps"],
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
        
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF \
        -dont-use-soft-clipped-bases \
        --native-pair-hmm-threads {threads} \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp {params.dbsnp} \
        {params.contigs} 

        """

rule combinevar:
    input:
        lambda wildcards: expand("analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz", sample=var_calling_units['sample'].unique(), contig_group=wildcards.contig_group)

    output:
        touch=touch("analysis/06_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/06_combinevar/{contig_group}.genomicsdb"),
    benchmark:
        "benchmarks/06_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs} 

        """

rule jointgeno:
    input:
        "analysis/06_combinevar/{contig_group}.done"
    output:
        vcf="analysis/07_jointgeno/all.{contig_group}.vcf.gz",
    benchmark:
        "benchmarks/07_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/06_combinevar/{contig_group}.genomicsdb"
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf} 
        """

rule merge_and_filter_vcf:
    input:
        expand("analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
    output:
        raw="analysis/08_merge_and_filter/all.merged.vcf.gz",
        filt="analysis/08_merge_and_filter/all.merged.filt.vcf.gz",
        pass_only="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz",
        vt_peek_raw="analysis/08_merge_and_filter/all.merged.vcf.gz.vt_peek.txt",
        vt_peek_pass="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz.vt_peek.txt"
    benchmark:
        "benchmarks/08_merge_and_filter/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        in_vcfs=expand("--INPUT analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
    envmodules:
        config['modules']['gatk'],
        config['modules']['vt']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --OUTPUT {output.raw} 

        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw}

        echo "mergeVcfs done."
        echo "mergeVcfs done." 1>&2

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {output.filt} 
        
        echo "VariantFiltration done." 
        echo "VariantFiltration done." 1>&2

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        -O {output.pass_only} 
        
        echo "SelectVariants done."
        echo "SelectVariants done." 1>&2

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 
        """

rule variant_annot:
    input:
        "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        html="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.html",
        vcf="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
        tbi="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
        html_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
        vcf_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        tbi_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
    benchmark:
        "benchmarks/09a_variant_annot/benchmark.txt"
    params:
        db_id=config["ref"]["snpeff_db_id"],
    envmodules:
        config['modules']['snpeff'],
        config['modules']['htslib']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards)
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -canon \
        -onlyProtein \
        -stats {output.html_canon} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf_canon}

        tabix {output.vcf_canon} 


        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -onlyProtein \
        -stats {output.html} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf}

        tabix {output.vcf} 
        """


rule snprelate:
    input:
        "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        "analysis/09b_snp_pca_and_dendro/report.html"
    params:
        gds="analysis/09b_snp_pca_and_dendro/all.gds"
    envmodules:
        config['modules']['R']
    conda:
        "envs/R.yaml"
    threads: 1
    resources:
        mem_gb = 60
    script:
        "bin/snprelate.Rmd"
