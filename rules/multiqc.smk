multiqc_input = []
if config["PE_or_SE"] =="SE":
    multiqc_input.append(expand("qc/fastqc/{units.sample}-SE_fastqc.html", units=units.itertuples()))
    multiqc_input.append(expand("qc/fastqc/{units.sample}-SE_fastqc.zip", units=units.itertuples()))
    multiqc_input.append(expand("trimmed_data/{units.sample}-SE.fastq.gz_trimming_report.txt", units=units.itertuples()))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))
elif config["PE_or_SE"] =="PE":
    multiqc_input.append(expand("qc/fastqc/{units.sample}-{read}_fastqc.html", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("qc/fastqc/{units.sample}-{read}_fastqc.zip", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("trimmed_data/{units.sample}-{read}.fastq.gz_trimming_report.txt", units=units.itertuples(), read=["R1","R2"]))
    multiqc_input.append(expand("analysis/star/{units.sample}.Log.final.out", units=units.itertuples()))

rule multiqc:
    input:
        multiqc_input
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
