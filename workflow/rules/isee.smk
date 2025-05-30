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

        """

rule deploy_isee_to_shinyappio:
    """
    Deploy iSEE to shinyapp.io and try to increase memory (may fail if account tier is too low).
    """
    input:
        sce="results/iSEE/sce.rds",
        app="results/iSEE/app.R",
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']))
    output:
        rsconnect=directory("results/iSEE/rsconnect/"), # meta data for this app for shinyapps.io
        deployed=touch("results/iSEE/deployed"),
        rscignore="results/iSEE/.rscignore"
    benchmark:
        "benchmarks/deploy_isee_to_shinyappio/bench.txt"
    params:
        deployed_file=lambda wildcards, output: os.path.basename(output.deployed),
        app_dir=lambda wildcards, input: f"appDir='{os.path.dirname(input.app)}'",
        app_info=f"appName = '{config['iSEE_app_name']}', account = '{config['shinyio_account_name']}'",
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock),
    envmodules:
        config['modules']['R'],
    threads: 1
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        echo '{params.deployed_file}' > {output.rscignore}

        Rscript --vanilla -e "renv::load('{params.renv_rproj_dir}'); options(repos = BiocManager::repositories()); rsconnect::deployApp({params.app_dir}, {params.app_info}); try(rsconnect::configureApp({params.app_dir}, {params.app_info}, size='xxxlarge'))"

        """
