#!/usr/bin/env Rscript
library(optparse)
suppressMessages(library(dplyr))
library(readr)
library(stringr)
 
option_list <- list(
  make_option(c("-s", "--sample_rgx"), type="character", default="^(\\S+)(?=_L000)", 
              help="regex for sample [default= %default]", metavar="character"),
  make_option(c("-r", "--group_rgx"), type="character", default="^(\\S+)(?=_L000)", 
              help="regex for group [default= %default]", metavar="character"),
  make_option(c("-g", "--geno_rgx"), type="character", default="^(\\S+)(?=_L000)", 
              help="regex for genotype [default= %default]", metavar="character"),
  make_option(c("-c", "--cond_rgx"), type="character", default="^(\\S+)(?=_L000)", 
              help="regex for condition [default= %default]", metavar="character"),
  make_option(c("-u", "--unit_rgx"), type="character", default="^(\\S+)(?=_L000)", 
              help="regex for unit [default= %default]", metavar="character")
); 
 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


fq_files <- list.files("../raw_data/", pattern = "*fastq.gz", recursive=TRUE)
R1_files <- grep("_R1_", fq_files, value = TRUE)

# sample    group    genotype    condition    unit    fq1    fq2    strandedness
df <- data.frame(fq1 = R1_files) %>%
mutate(sample = str_extract(basename(fq1), opt$sample_rgx),
       group = str_extract(basename(fq1), opt$group_rgx),
       genotype = str_extract(basename(fq1), opt$geno_rgx),
       condition = str_extract(basename(fq1), opt$cond_rgx),
       unit = str_extract(basename(fq1), opt$unit_rgx),
       fq2 = str_replace(fq1, "_R1_", "_R2_"),
       strandedness = "reverse") %>%
select(sample,group,genotype,condition,unit,fq1,fq2,strandedness)

# make sure no fq file listed more than once.
stopifnot(length(df$fq1) == length(unique(df$fq1)))
stopifnot(length(df$fq2) == length(unique(df$fq2)))

write_tsv(df, "units_template.tsv")
