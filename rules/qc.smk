
rule fastqc_raw:
    input:
        "raw_reads/{sample}.fastq.gz"
    output:
        html="qc/fastqc/raw_reads/{sample}_fastqc.html",
        zip="qc/fastqc/raw_reads/{sample}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/raw_reads/{sample}.log"
    wrapper:
        "0.35.1/bio/fastqc"

rule multiqc:
    input:
        # fastqc raw data
        expand("qc/fastqc/raw_reads/{samples.sample}_{read}_fastqc.html", read=["1","2"], samples=samples.itertuples()),
        expand("qc/fastqc/raw_reads/{samples.sample}_{read}_fastqc.zip", read=["1","2"], samples=samples.itertuples()),
        # trim_galore data----
            # cutadapt report from trim_galore
        expand("trimmed_data/{samples.sample}_{read}.fastq.gz_trimming_report.txt", read=["1","2"], samples=samples.itertuples()),
            # fastqc on trim_galore
        # STAR alingments
        expand("star/{samples.sample}.Log.final.out", samples=samples.itertuples()),
    params:
        # skip the pass1 from STAR
        "--ignore '*._STARpass1/*'"
    output:
        "qc/multiqc_report.html",
        "qc/multiqc_report_data/multiqc.log",
        "qc/multiqc_report_data/multiqc_cutadapt.txt",
        "qc/multiqc_report_data/multiqc_fastqc.txt",
        "qc/multiqc_report_data/multiqc_general_stats.txt",
        "qc/multiqc_report_data/multiqc_sources.txt",
        "qc/multiqc_report_data/multiqc_star.txt",
    log:
        "logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"
