renv::load(snakemake@params[["renv_rproj_dir"]])
.libPaths()

sce_file <- snakemake@input[['sce']]
comps <- snakemake@params[['comparisons']]
de_res_outfiles_dir <- snakemake@params[['de_res_outfiles_dir']]
out_sce <- snakemake@output[['sce']]

options(warn=1) # show all warnings

library(dplyr)
library(SingleCellExperiment)
library(DESeq2)
library(readr)

sce <- readRDS(sce_file)

for (i in seq_along(comps)){
    curr_comp <- comps[i]
    message("Parsing DE results for ", curr_comp, ".")
    de_res <- readRDS(file.path(de_res_outfiles_dir, curr_comp, "de_res.rds"))
    rowData(sce)[[paste0(curr_comp, ".padj")]] <- de_res$padj[match(rowData(sce)$ens_gene, de_res$ens_gene)]
    rowData(sce)[[paste0(curr_comp, ".LFC")]] <- de_res$log2FoldChange[match(rowData(sce)$ens_gene, de_res$ens_gene)]
    rowData(sce)[[paste0(curr_comp, ".rank")]] <- rank(rowData(sce)[[paste0(curr_comp, ".padj")]], ties.method = "min")
}

write_rds(sce, out_sce)
