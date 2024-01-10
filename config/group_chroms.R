library(yaml)
library(readr)
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))

# params from workflow
organism <- read_yaml("config.yaml")$organism

if (organism == "human"){
	ref_fasta <- read_yaml("config.yaml")$ref$human$sequence
} else if (organism == "mouse"){
	ref_fasta <- read_yaml("config.yaml")$ref$mouse$sequence
} else if (organism == "fly"){
	ref_fasta <- read_yaml("config.yaml")$ref$fly$sequence
} else if (organism == "worm"){
	ref_fasta <- read_yaml("config.yaml")$ref$worm$sequence
} else if (organism == "other"){
	ref_fasta <- read_yaml("config.yaml")$ref$other$sequence
}

outfile <- "grouped_contigs.tsv"

# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")

# output the standard chromosome names as a space-separated text file
std_chroms <- standardChromosomes(ref_gr)
nonstd <- seqlevels(ref_gr)[which(!seqlevels(ref_gr) %in% std_chroms)]

std_df <- data.frame(name=std_chroms, contigs=std_chroms)
nonstd_df <- data.frame(name="unplaced_contigs", contigs=paste(nonstd, collapse=','))

rbind(std_df, nonstd_df) %>% write_tsv(outfile)
