rule edger:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam.bai", units=units.itertuples()),
        "qc/multiqc_report.html" # require multiQC to be run before this analysis
    output:
        "rules/diffExp.html",
        "deliverables/GeneCounts.tsv"
    singularity:
        "docker://dpettinga/bbcrna:latest"
    script:
        "diffExp.Rmd"
