rule mergeTECounts:
    input:
        expand("analysis/TEcount/{units.sample}.cntTable", units=units.itertuples()),
    output:
        "deliverables/GeneTEcounts.tsv"
    log:
        "logs/mergeCounts.log"
    conda:
        "../envs/R.yaml"
    script:
        "mergeTEcounts.R"
