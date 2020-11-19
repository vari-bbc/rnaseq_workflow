import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.14.0")

##### check if conda is available. If not, we don't run rules requiring conda envs. Applies to dryruns and actual runs.
conda_avail = (which('conda') is not None)

##### load config and sample sheets #####

configfile: "bin/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

# units = pd.read_table("bin/units_test.tsv").set_index("sample", drop=False)
units = pd.read_table(config["units"]).set_index("sample", drop=False)
var_calling_units = pd.read_table("bin/variant_calling_units.tsv").set_index("unit", drop=False)

contrasts = pd.read_table(config["contrasts"]).set_index("name", drop=False)
#validate(contrasts, schema="schemas/contrasts.schema.yaml")

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


include: rules/variant_calling.smk

##### target rules #####

rule all:
    input:
        #lambda wildcards: expand("analysis/05_haplotypecaller/{units.sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.vcf.gz", units=var_calling_units.itertuples()) if config["call_variants"] else [],
        #lambda wildcards: expand("analysis/07_jointgeno/all.{contig_group}.vcf.gz", contig_group=contig_groups.name) if config["call_variants"] else [],
        #lambda wildcards: "analysis/08_merge_and_filter/all.merged.filt.vcf.gz" if config["call_variants"] else [],
        lambda wildcards: "analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi" if config["call_variants"] else [],
        "analysis/09b_snp_pca_and_dendro/report.html" if (config["call_variants"] and conda_avail) else [],
        # mergeLanesAndRename
            # SE
        # expand("raw_data/{units.sample}-SE.fastq.gz", units=units.itertuples()),
            # PE
        # expand("raw_data/{units.sample}-R1.fastq.gz", units=units.itertuples()),
        # expand("raw_data/{units.sample}-R2.fastq.gz", units=units.itertuples()),
        # fastq_screen
            # PE
        # expand("analysis/fastq_screen/{units.sample}-R1_screen.html", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R1_screen.txt", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R2_screen.html", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-R2_screen.txt", units=units.itertuples()),
            # SE
        # expand("analysis/fastq_screen/{units.sample}-SE_screen.html", units=units.itertuples()),
        # expand("analysis/fastq_screen/{units.sample}-SE_screen.txt", units=units.itertuples()),
        # Trim_Galore
            # SE
        # expand("analysis/trimmed_data/{units.sample}-SE_trimmed.fq.gz", units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-SE_fastqc.html", units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-SE_fastqc.zip", units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-SE.fastq.gz_trimming_report.txt", units=units.itertuples()),
            # PE
        # expand("analysis/trimmed_data/{units.sample}_R{read}_val_{read}.fq.gz", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.html", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.zip", read=[1,2], units=units.itertuples()),
        # expand("analysis/trimmed_data/{units.sample}-R{read}.fastq.gz_trimming_report.txt", read=[1,2], units=units.itertuples()),
        # STAR alignment
        # expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        # expand("analysis/star/{units.sample}.Log.out", units=units.itertuples()),
        # multiQC
        "analysis/multiqc/multiqc_report.html",
        #expand("analysis/02_splitncigar/{units.sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam", units=var_calling_units.itertuples())
        # edgeR
        #"bin/diffExp.html",

rule mergeLanesAndRename_SE:
    input:
    output:      "raw_data/{sample}-SE.fastq.gz"
    log:         "logs/mergeLanesAndRename/mergeLanesAndRename_SE-{sample}.log"
                 "logs/mergeLanesAndRename/mergeLanesAndRename_PE-{sample}.log"
    threads: 1
    resources:
        nodes =   1,
        mem_gb =  16,
    envmodules:  "bbc/R/R-3.6.0"
    script:      "bin/mergeLanesAndRename.R"

rule mergeLanesAndRename_PE:
    input:
    output:      "raw_data/{sample}-R1.fastq.gz",
                 "raw_data/{sample}-R2.fastq.gz"
    log:
                 "logs/mergeLanesAndRename/mergeLanesAndRename_PE-{sample}.log"
    threads: 1
    resources:
        nodes =   1,
        mem_gb =  16,
    envmodules:  "bbc/R/R-3.6.0"
    script:      "bin/mergeLanesAndRename.R"

def fastq_screen_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "raw_data/{sample}-SE.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 =    "raw_data/{sample}-R1.fastq.gz".format(**wildcards)
        R2 =    "raw_data/{sample}-R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule fastq_screen_PE:
    input:
        R1 =      "raw_data/{sample}-R1.fastq.gz",
        R2 =      "raw_data/{sample}-R2.fastq.gz",
    output:
        R1_html = "analysis/fastq_screen/{sample}-R1_screen.html",
        R1_txt =  "analysis/fastq_screen/{sample}-R1_screen.txt",
        R2_html = "analysis/fastq_screen/{sample}-R2_screen.html",
        R2_txt =  "analysis/fastq_screen/{sample}-R2_screen.txt",
    log:
        R1 =      "logs/fastq_screen/fastq_screen.{sample}-R1.log",
        R2 =      "logs/fastq_screen/fastq_screen.{sample}-R2.log",
    benchmark:    "benchmarks/fastq_screen/{sample}.bmk"
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  64,
    envmodules:   "bbc/fastq_screen/fastq_screen-0.14.0"
    shell:
        """
        fastq_screen --outdir analysis/fastq_screen/ {input.R1} 2> {log.R1}
        fastq_screen --outdir analysis/fastq_screen/ {input.R2} 2> {log.R2}
        """

rule fastq_screen_SE:
    input:
                    "raw_data/{sample}-SE.fastq.gz",
    output:
        html =      "analysis/fastq_screen/{sample}-SE_screen.html",
        txt =       "analysis/fastq_screen/{sample}-SE_screen.txt",
    log:
                    "logs/fastq_screen/fastq_screen.{sample}-SE.log",
    benchmark:
                    "benchmarks/fastq_screen/fastq_screen.{sample}.bmk"
    threads: 8
    resources:
        nodes =     1,
        mem_gb =    64,
    envmodules:     "bbc/fastq_screen/fastq_screen-0.14.0"
    shell:
        """
        fastq_screen --outdir analysis/fastq_screen/ {input} 2> {log.R1}
        """

def trim_galore_input(wildcards):
    if config["PE_or_SE"] == "SE":
        reads = "raw_data/{sample}-SE.fastq.gz".format(**wildcards)
        return reads
    elif config["PE_or_SE"] == "PE":
        R1 = "raw_data/{sample}-R1.fastq.gz".format(**wildcards)
        R2 = "raw_data/{sample}-R2.fastq.gz".format(**wildcards)
        return [R1,R2]

rule trim_galore_PE:
    input:
        trim_galore_input
    output:
        "analysis/trimmed_data/{sample}-R1_val_1.fq.gz",
        "analysis/trimmed_data/{sample}-R1_val_1_fastqc.html",
        "analysis/trimmed_data/{sample}-R1_val_1_fastqc.zip",
        "analysis/trimmed_data/{sample}-R1.fastq.gz_trimming_report.txt",
        "analysis/trimmed_data/{sample}-R2_val_2.fq.gz",
        "analysis/trimmed_data/{sample}-R2_val_2_fastqc.html",
        "analysis/trimmed_data/{sample}-R2_val_2_fastqc.zip",
        "analysis/trimmed_data/{sample}-R2.fastq.gz_trimming_report.txt"
    params:
        outdir="analysis/trimmed_data/"
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        "bbc/trim_galore/trim_galore-0.6.0"
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  80,
    shell:
        """
        trim_galore \
        --paired \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q 20 \
        --fastqc \
        2> {log.stderr} 1> {log.stdout}
        """

rule trim_galore_SE:
    input:
        trim_galore_input
    output:
        "analysis/trimmed_data/{sample}-SE_trimmed.fq.gz",
        "analysis/trimmed_data/{sample}-SE_trimmed_fastqc.zip",
        "analysis/trimmed_data/{sample}-SE_trimmed_fastqc.html",
        "analysis/trimmed_data/{sample}-SE.fastq.gz_trimming_report.txt",
    params:
        outdir="analysis/trimmed_data/"
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        "bbc/trim_galore/trim_galore-0.6.0"
    threads: 4
    resources:
        nodes =   1,
        mem_gb =  80,
    shell:
        """
        trim_galore \
        {input} \
        --output_dir {params.outdir} \
        --cores {threads} \
        -q 20 \
        --fastqc \
        2> {log.stderr} 1> {log.stdout}
        """

def STAR_input(wildcards):
    if config["PE_or_SE"] == "SE":
        fq1="analysis/trimmed_data/{sample}-SE_trimmed.fq.gz".format(**wildcards)
        return fq1
    elif config["PE_or_SE"] == "PE":
        fq1 = "analysis/trimmed_data/{sample}-R1_val_1.fq.gz".format(**wildcards)
        fq2 = "analysis/trimmed_data/{sample}-R2_val_2.fq.gz".format(**wildcards)
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
    log:
        "logs/star/{sample}.log"
    benchmark:
        "benchmarks/star/{sample}.txt"
    envmodules:
        "bbc/STAR/STAR-2.7.3a",
        "bbc/samtools/samtools-1.9"
    threads: 8
    resources:
        nodes =   1,
        mem_gb =  120,
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
        --outStd Log 2> {log}

        samtools index {output.bam}
        """

multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("analysis/fastq_screen/{units.sample}-SE_screen.txt", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed.fq.gz", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed_fastqc.zip", units=units.itertuples()))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-SE_trimmed_fastqc.html", units=units.itertuples()))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("analysis/fastq_screen/{units.sample}-R{read}_screen.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}.fq.gz", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.html", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}_val_{read}_fastqc.zip", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/trimmed_data/{units.sample}-R{read}.fastq.gz_trimming_report.txt", units=units.itertuples(), read=["1","2"]))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))

rule multiqc:
    input:
        multiqc_input
    params:
        "analysis/star/",
        "analysis/trimmed_data/",
        "analysis/fastq_screen/",
    output:
        "analysis/multiqc/multiqc_report.html",
        "analysis/multiqc/multiqc_report_data/multiqc.log",
        "analysis/multiqc/multiqc_report_data/multiqc_cutadapt.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_fastqc.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_general_stats.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_sources.txt",
        "analysis/multiqc/multiqc_report_data/multiqc_star.txt",
    log:
        "logs/multiqc.log"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 32,
    envmodules:
        "bbc/multiqc/multiqc-1.8"
    shell:
        """
        multiqc -f {params} \
        -o analysis/multiqc \
        --ignore '*._STARpass1/*' \
        -n multiqc_report.html \
        --cl-config 'max_table_rows: 999999' \
        2> {log}
        """

rule edgeR:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        "analysis/multiqc/multiqc_report.html" # require multiQC to be run before this analysis
    output:
        "bin/diffExp.html"
    log:
        "logs/edgeR.log"
    benchmark:
        "benchmarks/edgeR/edgeR.txt"
    envmodules:
        #use node095 RStudio Server R install
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
    script:
        "bin/diffExp.Rmd"


