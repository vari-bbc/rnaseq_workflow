# logging--------------------------------------------------------------
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="message")

# load libraries-------------------------------------------------------
library(vroom)
library(tidyverse)
library(biomaRt)

# load biomaRt---------------------------------------------------------
mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                         host="asia.ensembl.org",
                         dataset = 'hsapiens_gene_ensembl',
                         ensemblRedirect = FALSE)
t2g <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"),mart=mart) %>%
  dplyr::rename(gene_or_TE=ensembl_gene_id)

# import data----------------------------------------------------------
files <- snakemake@input
# initiate empty tibble with just gene/TE col
counts <- vroom(snakemake@input[[1]]) %>% dplyr::select(`gene/TE`)
# add all files to 'counts'
for (file in snakemake@input){
  counts <- vroom(file) %>% right_join(., counts)
}

# clean the data-------------------------------------------------------
# rename the colnames
colnames(counts) <- gsub("analysis/star/","",colnames(counts)) %>% 
  gsub(".Aligned.out.bam","",.) %>%
  gsub("gene/TE","gene_or_TE",.)
  
counts <- counts %>%
  mutate_at(vars(gene_or_TE), ~sub("(\\.(.+$))", "",.)) %>% # remove decimals from ENSG_IDs
  dplyr::left_join(.,t2g, by= "gene_or_TE") %>%
  dplyr::select(c("gene_or_TE","external_gene_name",everything()))

# write the output
write_tsv(counts,snakemake@output[[1]])

# save .RData allows debugging with snakemake object
save.image(file="rules/mergeCounts.RData")