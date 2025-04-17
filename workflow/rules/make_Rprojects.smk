rule make_Rproject:
    input:
        R_pkges="config/R_proj_packages.txt"
    output:
        rproj = "results/{Rproj}/{Rproj}.Rproj",
        renv_lock = "results/{Rproj}/renv.lock",
        renv_dependencies = "results/{Rproj}/_dependencies.R",
        renviron = "results/{Rproj}/.Renviron",
        rprofile = "results/{Rproj}/.Rprofile",
        renv_dir = directory("results/{Rproj}/renv"),
    benchmark:
        "benchmarks/make_Rproject/{Rproj}.txt"
    params:
        Rproj_dir = lambda wildcards,output: os.path.dirname(output.rproj),
        Renv_root = config['renv_root_path'],
        # tell renv not to use cache (actually install packages in project instead of symlinking) if requested in config file.
        set_use_renv_cache = lambda wildcards, output: "" if config['renv_use_cache'] else "printf \"RENV_CONFIG_CACHE_ENABLED=FALSE\\n\" >> {renviron}".format(renviron=output.renviron),
        # tell renv not to copy from user library if requested in config file.
        set_use_user_lib = lambda wildcards, output: "" if config['renv_use_user_lib'] else "printf \"RENV_CONFIG_INSTALL_SHORTCUTS=FALSE\\n\" >> {renviron}".format(renviron=output.renviron),
        set_symlink_from_cache = lambda wildcards, output: "" if config['renv_symlink_from_cache'] else "printf \"RENV_CONFIG_CACHE_SYMLINKS=FALSE\\n\" >> {renviron}".format(renviron=output.renviron),
        snakemake_dir=snakemake_dir,
    threads: 1
    envmodules:
        config['modules']['R']
    resources:
        mem_gb=64,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log", 
    shell:
        """
        Rscript --vanilla workflow/scripts/make_Rproject.R {params.Rproj_dir}
        
        # Make R script with package dependencies
        cat {input.R_pkges} | perl -lne "print qq:library(\$_):" > {output.renv_dependencies}
        
        # set RENV_PATHS_CACHE in .Renviron
        printf "RENV_PATHS_ROOT='{params.Renv_root}'\\n" > {output.renviron}

        {params.set_use_renv_cache}
        {params.set_use_user_lib}
        {params.set_symlink_from_cache}

        cd {params.Rproj_dir}
        
        # Bioc packages not already installed in user library can fail to install. Here we install problematic packages
        # See https://github.com/rstudio/renv/issues/81#issuecomment-497131224
        # Don't use --vanilla because we need to get the RENV_PATHS_ROOT from the .renviron file.
        Rscript {params.snakemake_dir}/workflow/scripts/install_renv_pkges.R {params.snakemake_dir}/{input.R_pkges}
        
        # Set up git
        git init
        git add -A 
        git commit -m "initial commit"
        
        """

