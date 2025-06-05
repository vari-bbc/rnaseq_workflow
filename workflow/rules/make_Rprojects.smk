rule make_Rproject:
    input:
        R_pkges="config/R_proj_packages.txt",
    output:
        rproj = "results/{Rproj}/{Rproj}.Rproj",
        gitignore = "results/{Rproj}/.gitignore",
        R_dir = directory("results/{Rproj}/R"),
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
        set_use_renv_cache = lambda wildcards, output: "RENV_CONFIG_CACHE_ENABLED={choice}".format(choice="TRUE" if config['renv_use_cache'] else "FALSE"),
        # tell renv not to copy from user library if requested in config file.
        set_use_user_lib = lambda wildcards, output: "RENV_CONFIG_INSTALL_SHORTCUTS={choice}".format(choice="TRUE" if config['renv_use_user_lib'] else "FALSE"),
        set_symlink_from_cache = lambda wildcards, output: "RENV_CONFIG_CACHE_SYMLINKS={choice}".format(choice="TRUE" if config['renv_symlink_from_cache'] else "FALSE"),
        set_use_pak = lambda wildcards, output: "RENV_CONFIG_PAK_ENABLED={choice}".format(choice="TRUE" if config['renv_use_pak'] else "FALSE"),
        install_pak = lambda wildcards, output: "library(pak)" if config['renv_use_pak'] else "",
        snakemake_dir=snakemake_dir,
        orgdb = config['orgdb'],
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
        echo '{params.install_pak}' > {output.renv_dependencies}
        (cat {input.R_pkges} && echo {params.orgdb}) | perl -lne 's:.*/::; s:@.*$::; print qq:library($_):' >> {output.renv_dependencies}
        
        # set RENV_PATHS_CACHE in .Renviron
        echo "RENV_PATHS_ROOT='{params.Renv_root}'" > {output.renviron}

        echo '{params.set_use_renv_cache}' >> {output.renviron}
        echo '{params.set_use_user_lib}' >> {output.renviron}
        echo '{params.set_symlink_from_cache}' >> {output.renviron}
        echo '{params.set_use_pak}' >> {output.renviron}


        cd {params.Rproj_dir}
       
        echo "Running script to set up renv library..."
        # Don't use --vanilla because we need to get the RENV_PATHS_ROOT from the .renviron file.
        Rscript {params.snakemake_dir}/workflow/scripts/install_renv_pkges.R {params.snakemake_dir}/{input.R_pkges}
        
        # Set up git
        echo "Setting up git now..."
        git init
        git add -A 
        git commit -m "initial commit"
        
        """

