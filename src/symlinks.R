library(readr)
library(dplyr)
# make symlinks
samples= read_delim("src/samples.tsv", delim="\t")

for (i in 1:nrow(samples)){
  if (!(paste0(samples[i,]$sample,"_1.fastq.gz") %in% list.files("raw_reads/"))){
    system(paste0("ln -s ",paste0(getwd(),"/raw_reads/",samples[i,]$fq1),
                  " ",
                  paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_1.fastq.gz")))
  }

  if (!(paste0(samples[i,]$sample,"_2.fastq.gz") %in% list.files("raw_reads/"))){
    system(paste0("ln -s ",paste0(getwd(),"/raw_reads/",samples[i,]$fq2),
                  " ",
                  paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_2.fastq.gz")))
  }
    # to UNLINK
  # system(paste0("unlink ",paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_1.fastq.gz")))
  # system(paste0("unlink ",paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_2.fastq.gz")))
  }
