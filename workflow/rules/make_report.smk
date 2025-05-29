rule make_final_report:
    input:
        website_template = expand("resources/report_template/{fname}", fname=["_site.yml", "footer.html", "index.Rmd", "multiqc.Rmd", "references.bib", "styles.css", "images/VAI_2_Line_White.png"]),
        de_res = expand("results/deseq2/DESeq2_{comp}.html", comp=pd.unique(comparisons["comparison_name"])),
        de_res_figs = expand("results/deseq2/deseq2_out_files/{comp}/individual_figures", comp=pd.unique(comparisons["comparison_name"])),
        de_res_tables = expand("results/deseq2/deseq2_out_files/{comp}/de_res.tsv", comp=pd.unique(comparisons["comparison_name"])),
        gsea = expand("results/gsea/gsea_{comp}.html", comp=pd.unique(comparisons["comparison_name"])),
        gsea_figs = expand("results/gsea/{comp}_out_files/individual_figures", comp=pd.unique(comparisons["comparison_name"])),
        gsea_tables = expand("results/gsea/{comp}_out_files/{collection}_gsea.xlsx", comp=pd.unique(comparisons["comparison_name"]), collection = config['pathway_str'].split(',')),
        renv_lock = ancient("results/{Rproj}/renv.lock".format(Rproj=config['Rproj_dirname'])),
        multiqc = "results/multiqc/multiqc_report.html"
    output:
        site_files = expand("results/make_final_report/{fname}", fname=["_site.yml", "footer.html", "index.Rmd","references.bib", "styles.css", "images/VAI_2_Line_White.png"]),
        multiqc_rmd = "results/make_final_report/multiqc.Rmd",
        de_res = expand("results/make_final_report/external_reports/DESeq2_{comp}.html", comp=pd.unique(comparisons["comparison_name"])),
        de_res_rmd = expand("results/make_final_report/DESeq2_{comp}.Rmd", comp=pd.unique(comparisons["comparison_name"])),
        de_res_figs = directory(expand("results/make_final_report/extras/deseq2_figures/{comp}", comp=pd.unique(comparisons["comparison_name"]))),
        de_res_tables = expand("results/make_final_report/extras/deseq2_tables/{comp}_de_res.tsv", comp=pd.unique(comparisons["comparison_name"])),
        gsea = expand("results/make_final_report/external_reports/gsea_{comp}.html", comp=pd.unique(comparisons["comparison_name"])),
        gsea_rmd = expand("results/make_final_report/gsea_{comp}.Rmd", comp=pd.unique(comparisons["comparison_name"])),
        gsea_figs = directory(expand("results/make_final_report/extras/gsea_figures/{comp}", comp=pd.unique(comparisons["comparison_name"]))),
        gsea_tables = expand("results/make_final_report/extras/gsea_tables/{comp}/{collection}_gsea.xlsx", comp=pd.unique(comparisons["comparison_name"]), collection = config['pathway_str'].split(',')),
        multiqc = "results/make_final_report/external_reports/multiqc_report.html",
        website = directory("results/make_final_report/BBC_RNAseq_Report"),
        isee_rmd = "results/make_final_report/iSEE.Rmd",
    benchmark:
        "benchmarks/make_final_report/bench.txt"
    params:
        template_dir = lambda wildcards, input: os.path.commonprefix(input.website_template),
        template_files = lambda wildcards, input: [ fname.replace("resources/report_template/", "")  for fname in input.website_template],
        de_res_comps = " ".join(["'" + comp + "'" for comp in pd.unique(comparisons["comparison_name"])]),
        de_res_yml = "\\n".join([f"        - text: {comp}\\n          href: DESeq2_{comp}.html" for comp in pd.unique(comparisons["comparison_name"])]),
        gsea_yml = "\\n".join([f"        - text: {comp}\\n          href: gsea_{comp}.html" for comp in pd.unique(comparisons["comparison_name"])]),
        ext_reports_dir = lambda wildcards, output: os.path.dirname(output.multiqc),
        renv_rproj_dir = lambda wildcards, input: os.path.dirname(input.renv_lock),
        root_dir = lambda wildcards, output: os.path.join(os.getcwd(), os.path.dirname(output.website)),
        de_res_outdir = lambda wildcards, input: os.path.dirname(os.path.dirname(input.de_res_figs[0])),
        gsea_dir = lambda wildcards, input: os.path.dirname(input.gsea[0]),
        isee_site_yml = "    - text: iSEE\\n      icon: fa-solid fa-gem\\n      href: isee.html" if config['deploy_to_shinyio'] else "",
        isee_app_name = config['iSEE_app_name'],
        shinyio_account = config['shinyio_account_name']
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
        cd {params.template_dir}
        cp --parents -r {params.template_files} {params.root_dir}/
        cd -

        perl -i -lnpe 's/<<<DE_RES>>>/{params.de_res_yml}/' {params.root_dir}/_site.yml
        perl -i -lnpe 's/<<<GSEA_RES>>>/{params.gsea_yml}/' {params.root_dir}/_site.yml
        perl -i -lnpe 's/<<<ISEE>>>/{params.isee_site_yml}/' {params.root_dir}/_site.yml

        ln -sr {input.multiqc} {output.multiqc}
        ln -sr {input.de_res} {params.ext_reports_dir}
        ln -sr {input.gsea} {params.ext_reports_dir}
       
        cat {output.multiqc_rmd} | perl -lnpe "s:MultiQC:iSEE:; s|./external_reports/multiqc_report.html|https://{params.shinyio_account}.shinyapps.io/{params.isee_app_name}|" > {output.isee_rmd}

        for comp in {params.de_res_comps}
        do
            mkdir -p {params.root_dir}/extras/deseq2_figures/${{comp}}/
            ln -sr {params.de_res_outdir}/${{comp}}/individual_figures/* {params.root_dir}/extras/deseq2_figures/${{comp}}/
            
            mkdir -p {params.root_dir}/extras/deseq2_tables/
            ln -sr {params.de_res_outdir}/${{comp}}/de_res.tsv {params.root_dir}/extras/deseq2_tables/${{comp}}_de_res.tsv

            cat {output.multiqc_rmd} | perl -lnpe "s:MultiQC:${{comp}}:; s:multiqc_report:DESeq2_${{comp}}:" > "{params.root_dir}/DESeq2_${{comp}}.Rmd"
            
            mkdir -p {params.root_dir}/extras/gsea_figures/${{comp}}/
            ln -sr {params.gsea_dir}/${{comp}}_out_files/individual_figures/* {params.root_dir}/extras/gsea_figures/${{comp}}/
            
            mkdir -p {params.root_dir}/extras/gsea_tables/
            ln -sr {params.gsea_dir}/${{comp}}_out_files/*_gsea.xlsx {params.root_dir}/extras/gsea_tables/${{comp}}/
            
            cat {output.multiqc_rmd} | perl -lnpe "s:MultiQC:${{comp}}:; s:multiqc_report:gsea_${{comp}}:" > "{params.root_dir}/gsea_${{comp}}.Rmd"
        done
        

        Rscript -e "renv::load('{params.renv_rproj_dir}'); rmarkdown::render_site('{params.root_dir}')" 
        """
