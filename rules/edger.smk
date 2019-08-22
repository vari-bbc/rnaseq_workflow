rule edgeR_longReport:
    input:
        "deliverables/UniquelyMappingRates.txt",
        "deliverables/UniquelyMappingReads.txt",
        "deliverables/starMatrix.txt",
        # sp = config["species"]["short"],
        # species = config["species"]["long"],
    output:
        # R Objects - also used for edgeR_shortReport
        expand("deliverables/r_objects/topGenes.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/r_objects/volcano.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/r_objects/contrastMatrix.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/r_objects/DGE_design.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        "deliverables/r_objects/pca_plot.rds",
        "deliverables/r_objects/var_plot.rds",
        "deliverables/r_objects/mat.rds",
        "deliverables/r_objects/rowNames_use.rds",
        "deliverables/r_objects/filtered_meta_heatmap.rds",
        # edgeR results
        expand("deliverables/{contrast}.txt", contrast=contrasts['contrast'].tolist()),
        directory("edgeR_longReport_cache/"),
        directory("edgeR_longReport_files/"),
        # GSEA results
        expand("deliverables/{contrast}.cls", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/{contrast}.gct", contrast=contrasts['contrast'].tolist()),
        # HTML report
        "deliverables/edgeR_longReport.html",
    conda:
        "../envs/edger.yaml"
    shell:
        '''
        Rscript -e 'rmarkdown::render("src/edgeR_longReport.Rmd",output_format="html_document")'
        mv src/edgeR_longReport.html deliverables/edgeR_longReport.html
        '''

rule edgeR_shortReport:
    input:
        # heatmap objects
        "deliverables/r_objects/pca_plot.rds",
        "deliverables/r_objects/var_plot.rds",
        "deliverables/r_objects/mat.rds",
        "deliverables/r_objects/rowNames_use.rds",
        "deliverables/r_objects/filtered_meta_heatmap.rds",
        # DE objects
        expand("deliverables/r_objects/DGE_design.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/r_objects/topGenes.{contrast}.rds", contrast=contrasts['contrast'].tolist()),
        expand("deliverables/r_objects/volcano.{contrast}.rds", contrast=contrasts['contrast'].tolist())
    output:
        "deliverables/edgeR_shortReport.html",
    # singularity:
    #     "shub://deanpettinga/rnaseq:edger_2"
    conda:
        "../envs/edger.yaml"
    shell:
        '''
        Rscript -e 'rmarkdown::render("src/edgeR_shortReport.Rmd",output_format="html_document")'
        mv src/edgeR_shortReport.html deliverables/edgeR_shortReport.html
        '''
