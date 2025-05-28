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

        """

rule deploy_isee_to_shinyappio:
    input:
        sce="results/iSEE/sce.rds",
        app="results/iSEE/app.R",
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']))
    output:
        deployed=touch("results/iSEE/deployed"),
        rscignore="results/iSEE/.rscignore"
    benchmark:
        "benchmarks/deploy_isee_to_shinyappio/bench.txt"
    params:
        deployed_file=lambda wildcards, output: os.path.basename(output.deployed),
        app_dir=lambda wildcards, input: os.path.dirname(input.app),
        R_proj_path = lambda wildcards, input: os.path.join(os.getcwd(), os.path.dirname(input.renv_lock)),
        isee_app_name = config['iSEE_app_name'],
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock),
        shinyio_account = config['shinyio_account_name'],
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

        Rscript --vanilla -e "renv::load('{params.renv_rproj_dir}'); options(repos = BiocManager::repositories()); rsconnect::deployApp(appDir='{params.app_dir}', appName = '{params.isee_app_name}', account = '{params.shinyio_account}')"

        """
