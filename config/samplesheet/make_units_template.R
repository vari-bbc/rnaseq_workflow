#!/usr/bin/env Rscript
library(optparse)
suppressMessages(library(dplyr))
library(readr)
library(stringr)
 
option_list <- list(
  make_option(c("-s", "--sample_rgx"), type="character", default="^([^_]+)", 
              help="regex for sample [default= %default]", metavar="character"),
  make_option(c("-r", "--group_rgx"), type="character", default="^([^_]+)", 
              help="regex for group [default= %default]", metavar="character")
); 
 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


fq_files <- list.files("../../raw_data/", pattern = "*fastq.gz", recursive=TRUE)
R1_files <- grep("_R1[_\\.]", fq_files, value = TRUE)

# sample    group    genotype    condition    unit    fq1    fq2    strandedness
df <- data.frame(fq1 = R1_files) %>%
mutate(sample = str_extract(basename(fq1), opt$sample_rgx),
       group = str_extract(basename(fq1), opt$group_rgx),
       fq2 = str_replace(fq1, "_R1([_\\.])", "_R2\\1"),
       RG="") %>%
arrange(sample) %>%
select(sample,group,fq1,fq2,RG)

# make sure no fq file listed more than once.
stopifnot(length(df$fq1) == length(unique(df$fq1)))
stopifnot(length(df$fq2) == length(unique(df$fq2)))
stopifnot(length(c(df$fq1, df$fq2)) == length(unique(c(df$fq1, df$fq2))))

# make sure all R2 files have 'R2' in the name
stopifnot(sum(str_detect(df$fq2, "_R2[_\\.]")) == length(df$fq2))
stopifnot(sum(str_detect(df$fq1, "_R1[_\\.]")) == length(df$fq1))

# make sure all found fastq files listed
stopifnot(all(sort(c(df$fq1, df$fq2)) == sort(fq_files)))

write_tsv(df, "units_template.tsv")
