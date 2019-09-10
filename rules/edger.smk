rule edgeR_longReport:
    input:
        "deliverables/UniquelyMappingRates.txt",
        "deliverables/UniquelyMappingReads.txt",
        "deliverables/starMatrix.txt",
        "deliverables/counts.tsv",
    output:
        # R Objects - also used for edgeR_shortReport
        expand("src/r_objects/topGenes.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("src/r_objects/volcano.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("src/r_objects/contrastMatrix.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("src/r_objects/DGE_design.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        "src/r_objects/pca_plot.rds",
        "src/r_objects/var_plot.rds",
        "src/r_objects/mat.rds",
        "src/r_objects/rowNames_use.rds",
        "src/r_objects/filtered_meta_heatmap.rds",
        # edgeR results
        expand("deliverables/{contrast}.txt", contrast=contrasts['contrast'].tolist()),
        "deliverables/cpm.txt",
        "deliverables/fpkm.txt",
        directory("src/edgeR_longReport_cache/"),
        directory("src/edgeR_longReport_files/"),
        # HTML report
        "edgeR_longReport.html",
    # conda:
    #     "../envs/edger.yaml"
    script:
        "src/edgeR_longReport.Rmd"
