library(yaml)
library(readr)
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))

# params from workflow
# check if quick_ref is specified
quick_ref <- read_yaml("config.yaml")$quick_ref
if (is.null(quick_ref$species_name)) {
    cat("Using the index files manually specified in the config file.\n")
    ref_fasta <- read_yaml("config.yaml")$ref$sequence
} else {
    cat("Using quick_ref instead of manually indicated index files.\n")
    if (is.null(quick_ref$ref_genome_version)) {
        cat("Version number not specified. Will use the latest version of the BBC-curated reference files.\n")
        quick_ref$ref_genome_version <- "latest"
    }
    ref_fasta <- paste0("/varidata/research/projects/bbc/versioned_references/", 
                        quick_ref$ref_genome_version, 
                        "/data/", 
                        quick_ref$species_name, 
                        "/sequence/",
                        quick_ref$species_name,
                        ".fa")
    if (!file.exists(ref_fasta)) {
        stop(paste("The quick_ref reference fasta file does not exist:", ref_fasta))
    }
}

outfile <- "grouped_contigs.tsv"

# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")

# output the standard chromosome names as a space-separated text file
std_chroms <- standardChromosomes(ref_gr)
nonstd <- seqlevels(ref_gr)[which(!seqlevels(ref_gr) %in% std_chroms)]

std_df <- data.frame(name=std_chroms, contigs=std_chroms)
nonstd_df <- data.frame(name="unplaced_contigs", contigs=paste(nonstd, collapse=','))

if(nonstd_df$contigs != ""){
    outdf <- rbind(std_df, nonstd_df)
} else{
    outdf <- std_df
}

write_tsv(outdf, outfile)
