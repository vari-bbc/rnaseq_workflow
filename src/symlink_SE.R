library(readr)
library(dplyr)
# make symlinks
units<- read_delim("src/units.tsv", delim="\t")

for (i in 1:nrow(units)){
  # symlink the files in column fq1 to "sample_unit_R1.fastq.gz
  system(paste0("ln -s ",paste0(getwd(),"/raw_reads/",units[i,]$fq1),
                  " ",
                  paste0("raw_reads/",units[i,]$sample,"-",units[i,]$unit, ".fastq.gz")))
  # enable snakemake to delete this symlink.
  system(paste0("touch -h ","raw_reads/",units[i,]$sample,"-",units[i,]$unit, ".fastq.gz"))

}
