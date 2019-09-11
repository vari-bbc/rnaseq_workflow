library(readr)
library(dplyr)
# make symlinks
units<- read_delim("src/units.tsv", delim="\t")

for (i in 1:nrow(units)){
  if (!(paste0(units[i,]$sample,"_",
               units[i,]$unit,
               "_R1.fastq.gz") %in% list.files("raw_reads/"))){
    system(paste0("ln -s ",paste0(getwd(),"/raw_reads/",units[i,]$fq1),
                  " ",
                  paste0("raw_reads/",units[i,]$sample,"_",units[i,]$unit, "_R1.fastq.gz")))
  }

  if (!(paste0(units[i,]$sample,"_",
               units[i,]$unit,
               "_R2.fastq.gz") %in% list.files("raw_reads/"))){
    system(paste0("ln -s ",paste0(getwd(),"/raw_reads/",units[i,]$fq2),
                  " ",
                  paste0("raw_reads/",units[i,]$sample,"_",units[i,]$unit, "_R2.fastq.gz")))
  }
    # to UNLINK
  # system(paste0("unlink ",paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_1.fastq.gz")))
  # system(paste0("unlink ",paste0(getwd(),"/raw_reads/",samples[i,]$sample,"_2.fastq.gz")))
}
