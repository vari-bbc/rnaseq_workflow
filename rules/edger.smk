rule edger:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples()),
        "qc/multiqc_report.html" # require multiQC to be run before this analysis
    output:
        "rules/diffExp.html"
    singularity:
        "docker://dpettinga/bbcrna"
    script:
        "diffExp.Rmd"
