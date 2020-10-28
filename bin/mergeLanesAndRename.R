log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(yaml)
config <- yaml.load_file("bin/config.yaml")
units <- read_tsv(config$units)
samp = snakemake@output[[1]] %>% gsub("raw_data/","",.) %>% gsub("(-.{2}).fastq.gz","",.)

# samp = "RNANVSA"
# units <- read_tsv("../bin/units.tsv")

# get 'units' for each sample
samp_units <- units %>% dplyr::filter(sample==samp)

# check that fastq files are gzipped (have .gz suffix)
fq1_ends_with_gz <- all(samp_units$fq1 %>% stringr::str_detect(".gz$"))
if(!fq1_ends_with_gz){
    stop(paste0("fq1 for ", samp, " not a gzipped file"))
}

# check fq2 similarly only if the config file indicates PE reads
if (config$PE_or_SE=="PE"){
    fq2_ends_with_gz <- all(samp_units$fq2 %>% stringr::str_detect(".gz$"))
    if(!fq2_ends_with_gz){
        stop(paste0("fq2 for ", samp, " not a gzipped file"))
    }
}


# Now merge lanes

# if this is a SE experiment
if (config$PE_or_SE=="SE"){
  print(paste("merging SE reads for:", samp))


  if (all(file.exists(paste0("raw_data/", samp_units$fq1)))){
    # merge SE reads
    bash_call <- ifelse(nrow(samp_units) > 1, 
                        paste0("cat ", gsub(toString(paste0("raw_data/", samp_units$fq1)), pattern=",", replacement = ""), " > "), 
                        paste0("ln -s ", samp_units$fq1, " "))
    bash_call <- paste0(bash_call, "raw_data/", samp, "-SE.fastq.gz")
    print(bash_call)
    system(bash_call)
    print("R1 units merged")
    save.image(file=paste0("logs/mergeLanesAndRename/mergeLanesAndRename_SE-",samp,".RData"))
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp, "listed in bin/units.tsv not present in raw_data/."))
  }

} else if (config$PE_or_SE=="PE"){
# if this is a PE experiment
  # merge R1 reads
  print(paste("merging PE reads for:", samp))


  if (all(file.exists(paste0("raw_data/", samp_units$fq1)))){
    bash_call <- ifelse(nrow(samp_units) > 1, 
                        paste0("cat ", gsub(toString(paste0("raw_data/", samp_units$fq1)), pattern=",", replacement = ""), " > "), 
                        paste0("ln -s ", samp_units$fq1, " "))
    bash_call <- paste0(bash_call, "raw_data/", samp, "-R1.fastq.gz")
    print(bash_call)
    system(bash_call)
    print("R1 units merged")
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq1 file for", samp, "listed in bin/units.tsv not present in raw_data/."))
  }
  # merge R2 reads
  if (all(file.exists(paste0("raw_data/", samp_units$fq2)))){
    bash_call <- ifelse(nrow(samp_units) > 1, 
                        paste0("cat ", gsub(toString(paste0("raw_data/", samp_units$fq2)), pattern=",", replacement = ""), " > "), 
                        paste0("ln -s ", samp_units$fq2, " "))
    bash_call <- paste0(bash_call, "raw_data/", samp, "-R2.fastq.gz")
    print(bash_call)
    system(bash_call)
    print("R2 units merged")
    save.image(file=paste0("logs/mergeLanesAndRename/mergeLanesAndRename_PE-",samp,".RData"))
  } else {
    stop(paste("Error in mergeLanesAndRename.R: fq2 file for", samp, "listed in bin/units.tsv not present in raw_data/."))
  }
  } else {
  stop("Error in mergeLanesAndRename.R: Neither SE nor PE specified.")
}
print("merging complete, exiting script.")
