rule make_final_report:
    input:
        website_template = expand("resources/report_template/{fname}", fname=glob_wildcards("resources/report_template/{fname}").fname),
        renv_lock = "results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']),
        multiqc = "results/multiqc/multiqc_report.html"
    output:
        root_dir = directory("results/make_final_report"),
        multiqc = "results/make_final_report/external_reports/multiqc_report.html",
        website = directory("results/make_final_report/BBC_RNAseq_Report")
    benchmark:
        "benchmarks/make_final_report/bench.txt"
    params:
        website_dir = lambda wildcards, output: os.path.dirname(output.website), 
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock)
    envmodules:
        config['modules']['R'],
        config['modules']['pandoc']
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        cp -r {input.website_template} {params.website_dir}
        cp {input.multiqc} {output.multiqc}
        
        
        Rscript -e "renv::load('{params.renv_rproj_dir}'); rmarkdown::render_site('{params.website_dir}')" 
        """
