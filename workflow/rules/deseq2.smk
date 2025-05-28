rule deseq2:
    input:
        se="results/SummarizedExperiment/SummarizedExperiment.rds",
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']))
    output:
        rmd="results/deseq2/DESeq2_{comparison}.Rmd",
        html="results/deseq2/DESeq2_{comparison}.html",
        de_res="results/deseq2/deseq2_out_files/{comparison}/de_res.tsv",
        de_res_rds="results/deseq2/deseq2_out_files/{comparison}/de_res.rds",
        counts="results/deseq2/deseq2_out_files/{comparison}/counts.tsv",
        vst="results/deseq2/deseq2_out_files/{comparison}/vst.tsv",
        vsd_rds="results/deseq2/deseq2_out_files/{comparison}/vsd.rds",
        fig_dir=directory("results/deseq2/deseq2_out_files/{comparison}/individual_figures"),
    benchmark:
        "benchmarks/deseq2/{comparison}.txt"
    params:
        deseq_template = "resources/deseq_template.Rmd",
        fdr_cutoff = config['fdr_cutoff'],
        comparison = lambda wildcards: wildcards.comparison,
        group_test = lambda wildcards: comparisons.loc[comparisons['comparison_name'] == wildcards.comparison, 'group_test'].values[0],
        group_reference = lambda wildcards: comparisons.loc[comparisons['comparison_name'] == wildcards.comparison, 'group_reference'].values[0],
        genes_of_interest = config['genes_of_interest'],
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
        cp {params.deseq_template} {output.rmd}

        Rscript --vanilla -e "renv::load('{params.renv_rproj_dir}'); rmarkdown::render('{output.rmd}', 
                              params = list(se_obj = '{input.se}',
                                            comparison_name = '{params.comparison}',
                                            group_test = '{params.group_test}', 
                                            group_reference = '{params.group_reference}',
                                            fdr_cutoff = {params.fdr_cutoff},
                                            genes_of_interest = '{params.genes_of_interest}'))"
        """

rule add_DE_to_SE:
    input:
        sce="results/SummarizedExperiment/sce.rds",
        res_rds=expand("results/deseq2/deseq2_out_files/{comparison}/de_res.rds", comparison = pd.unique(comparisons["comparison_name"])),
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname'])),
    output:
        sce="results/add_DE_to_SE/sce.rds"
    benchmark:
        "benchmarks/add_DE_to_SE/bench.txt"
    params:
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock),
        de_res_outfiles_dir = lambda wildcards, input: os.path.dirname(os.path.dirname(input.res_rds[0])),
        comparisons = pd.unique(comparisons["comparison_name"]),
    envmodules:
        config['modules']['R'],
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    script:
        "../scripts/add_DE_to_SCE.R"
