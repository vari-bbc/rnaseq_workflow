rule edger:
    input:
        expand("analysis/star/{units.sample}.Aligned.sortedByCoord.out.bam", units=units.itertuples())
    output:
        "rules/diffExp.html"
    singularity:
        "docker://dpettinga/bbcrna"
    script:
        "diffExp.Rmd"
