rule fastqc:
    input:
        "raw_reads/{sample}_{unit}_R{read}.fastq.gz",
    output:
        html="qc/fastqc/{sample}_{unit}_R{read}_fastqc.html",
        zip="qc/fastqc/{sample}_{unit}_R{read}_fastqc.zip",
    params: ""
    log:
        "logs/fastqc/{sample}_{unit}_R{read}.log"
    wrapper:
        "0.35.1/bio/fastqc"

rule multiqc:
    input:
        # fastqc raw data
        expand("qc/fastqc/{units.sample}_{units.unit}_R{read}_fastqc.html", read=["1","2"], units=units.itertuples()),
        expand("qc/fastqc/{units.sample}_{units.unit}_R{read}_fastqc.zip", read=["1","2"], units=units.itertuples()),
        # trim_galore data----
            # cutadapt report from trim_galore
        expand("trimmed_data/{units.sample}_{units.unit}_R{read}.fastq.gz_trimming_report.txt", read=["1","2"], units=units.itertuples()),
            # fastqc on trim_galore
        # STAR alingments
        expand("analysis/star/{units.sample}_{units.unit}.Log.final.out", units=units.itertuples()),
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
