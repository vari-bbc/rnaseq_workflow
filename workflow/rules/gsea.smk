rule gsea:
    input:
        res_rds="results/deseq2/deseq2_out_files/{comparison}/de_res.rds",
        vsd_rds="results/deseq2/deseq2_out_files/{comparison}/vsd.rds",
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname']))
    output:
        rmd="results/gsea/gsea_{comparison}.Rmd",
        html="results/gsea/gsea_{comparison}.html",
        figs=directory("results/gsea/{comparison}_out_files/individual_figures"),
        outfiles=expand("results/gsea/{{comparison}}_out_files/{collection}_gsea.xlsx", collection = config['pathway_str'].split(','))
    benchmark:
        "benchmarks/gsea/{comparison}.txt"
    params:
        orgdb = config['orgdb'],
        kegg_org = config['kegg_org'],
        reactome_org = config['reactome_org'],
        msigdb_organism = config['msigdb_organism'],
        pathway_str = config['pathway_str'],
        gsea_template = "resources/gsea_template.Rmd",
        fdr_cutoff = config['fdr_cutoff'],
        comparison = lambda wildcards: wildcards.comparison,
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
        cp {params.gsea_template} {output.rmd}

        Rscript --vanilla -e "renv::load('{params.renv_rproj_dir}'); rmarkdown::render('{output.rmd}', params = list(orgdb = '{params.orgdb}', kegg_org = '{params.kegg_org}', comparison_name = '{params.comparison}', reactome_org = '{params.reactome_org}', msigdb_organism = '{params.msigdb_organism}', fdr_cutoff = '{params.fdr_cutoff}', de_res = '{input.res_rds}', vsd_rds = '{input.vsd_rds}', pathway_str = '{params.pathway_str}'))"
        """
