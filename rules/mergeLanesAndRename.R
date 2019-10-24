library(tidyverse)
library(yaml)

config <- yaml.load_file("src/config.yaml")
units <- read_delim(config$units,delim = "\t")
samp = snakemake@output[[1]] %>% gsub("raw_reads/","",.) %>% gsub("(-.{2}).fastq.gz","",.)
            
# Now merge lanes
print(paste("merging reads for:",samp))

# if this is a SE experiment
if (config$PE_or_SE=="SE"){
  # merge SE reads
  system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-SE.fastq.gz"))
} else if (config$PE_or_SE=="PE"){
# if this is a PE experiment
  # merge R1 reads
  system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-R1.fastq.gz"))
  # merge R2 reads
  system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq2)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-R2.fastq.gz"))
} else {
  stop("Error in mergeLanesAndRename.R. Neither SE nor PE specified.")
}
