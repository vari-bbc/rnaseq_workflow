log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(yaml)
config <- yaml.load_file("src/config.yaml")
units <- read_tsv(config$units)
samp = snakemake@output[[1]] %>% gsub("raw_reads/","",.) %>% gsub("(-.{2}).fastq.gz","",.)

# samp = "RNANVSA"
# units <- read_tsv("../src/units.tsv")

# check that fastq files are gzipped (have .gz suffix)
fq1_ends_with_gz <- (units %>% dplyr::filter(sample==samp))$fq1 %>% stringr::str_detect(".gz$")
if(!fq1_ends_with_gz){
    stop(paste0("fq1 for ", samp, " not a gzipped file"))
}

# check fq2 similarly only if the config file indicates PE reads
if (config$PE_or_SE=="PE"){
    fq2_ends_with_gz <- (units %>% dplyr::filter(sample==samp))$fq2 %>% stringr::str_detect(".gz$")
    if(!fq2_ends_with_gz){
        stop(paste0("fq2 for ", samp, " not a gzipped file"))
    }
}


# Now merge lanes

# if this is a SE experiment
if (config$PE_or_SE=="SE"){
  print(paste("merging SE reads for:", samp))
  if (file.exists(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)))){
    # merge SE reads
    system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-SE.fastq.gz"))
    print("R1 units merged")
    save.image(file=paste0("logs/mergeLanesAndRename/mergeLanesAndRename_SE-",samp,".RData"))
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp, "listed in src/units.tsv not present in raw_reads/."))
  }
  
} else if (config$PE_or_SE=="PE"){
# if this is a PE experiment
  # merge R1 reads
  print(paste("merging PE reads for:", samp))
  if (file.exists(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)))){
    system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq1)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-R1.fastq.gz"))
    print("R1 units merged")
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp, "listed in src/units.tsv not present in raw_reads/."))
  }
  # merge R2 reads
  if (file.exists(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq2)))){
    system(paste0("cat ",
                gsub(toString(paste0("raw_reads/",(units %>% dplyr::filter(sample==samp))$fq2)),pattern=",", replacement = ""),
                " > ", "raw_reads/", samp, "-R2.fastq.gz"))
    print("R2 units merged")
    save.image(file=paste0("logs/mergeLanesAndRename/mergeLanesAndRename_PE-",samp,".RData"))
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq2 file for", samp, "listed in src/units.tsv not present in raw_reads/."))
  } 
  } else {
  stop("Error in mergeLanesAndRename.R: Neither SE nor PE specified.")
}
print("merging complete, exiting script.")
