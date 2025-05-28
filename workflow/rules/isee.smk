rule isee:
    input:
        sce="results/add_DE_to_SE/sce.rds",
        app="resources/iSEE_app.R",
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']))
    output:
        sce="results/iSEE/sce.rds",
        app="results/iSEE/app.R",
    benchmark:
        "benchmarks/isee/isee.txt"
    params:
        R_proj_path = lambda wildcards, input: os.path.join(os.getcwd(), os.path.dirname(input.renv_lock)),
        root_dir = lambda wildcards, output: os.path.join(os.getcwd(), os.path.dirname(output.app)),
        isee_app_name = config['iSEE_app_name']
    envmodules:
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        cp {input.sce} {output.sce}
        cp {input.app} {output.app}

        # cat {input.app} | perl -npe "s:<<<R_proj_path>>>:'{params.R_proj_path}':; s:<<<shiny_app_dir>>>:'{params.root_dir}':; s:<<<isee_app_name>>>:'{params.isee_app_name}':" > {output.app}
               
        """
