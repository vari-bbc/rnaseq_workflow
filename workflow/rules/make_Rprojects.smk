rule make_Rproject:
    input:
        se="results/SummarizedExperiment/SummarizedExperiment.rds",
    output:
        rproj= "{Rproj}/{Rproj}.Rproj"
    benchmark:
        "benchmarks/make_Rproject/{Rproj}.txt"
    params:
        Rproj_dir = lambda wildcards,output: os.path.dirname(output.rproj),
        R_packages = Rproj_packages
    threads: 1
    envmodules:
        config['modules']['R']
    resources:
        mem_gb=64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        Rscript --vanilla workflow/scripts/make_Rproject.R {params.Rproj_dir}
        
        cd {params.Rproj_dir}

        git init

        Rscript --vanilla -e "BiocManager::install({params.R_packages}); renv::snapshot();"

        # symlink unix results to R project folder
        mkdir unix_res
        cd unix_res
        ln -sr ../../{input.se} .
    """

