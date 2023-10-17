args <- commandArgs(trailingOnly = TRUE)


gtf_file <- args[1]
orgdb <- args[2]
out_se <- args[3]
out_sce <- args[4]
out_sizeFactors <- args[5]

star_dir <- "results/star"
salmon_dir <- "results/salmon"
samplesheet <- "config/samplesheet/units.tsv"

# Packages loaded

library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(edgeR)
# load the org.db for your organism
if(!require(orgdb, character.only=TRUE)){
    BiocManager::install(orgdb)
    library(orgdb, character.only=TRUE)
}
library(DESeq2)
library(AnnotationDbi)
library(tximport)
library(GenomicFeatures)
library(scater)
library(rjson)

# Infer the strand from Salmon results
# ISR or SR = reverse
# ISF or SF = forward
# IU or U = unstranded
salmon_lib_files <- list.files(salmon_dir, pattern = "lib_format_counts.json", recursive=TRUE, full.names = TRUE)
lib_types <- unlist(lapply(salmon_lib_files, function(x) fromJSON(file=x)$expected_format))
lib_type <- unique(lib_types)

if(length(lib_type) > 1){
 stop("More than 1 library type detected.")
}

if(!lib_type %in% c("ISR", "SR", "ISF", "SF", "IU", "U")){
 stop("Unknown library type detected.")
}

# Read counts
files <- list.files(star_dir, pattern = "ReadsPerGene.out.tab", full.names = FALSE)
names(files) <- str_remove_all(files, ".ReadsPerGene.out.tab$")

counts_col <- case_when(
    lib_type %in% c("ISR", "SR") ~ 4,
    lib_type %in% c("ISF", "SF") ~ 3,
    lib_type %in% c("IU", "U") ~ 2,
    .default=999
)
stopifnot(counts_col %in% c(2, 3, 4))

dge <- edgeR::readDGE(files, path = star_dir, columns = c(1, counts_col), 
                      skip=4, labels = names(files), header=FALSE)


raw_cts <- edgeR::getCounts(dge)
names(attributes(raw_cts)$dimnames) <- NULL

# Read TPMs

txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

files <- list.files(salmon_dir, recursive=TRUE, pattern = "quant.sf", full.names = TRUE)
names(files) <- basename(str_remove(files, "\\/quant.sf"))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
tpms <- txi.salmon$abundance

# some genes are only in the STAR counts
tpms <- tpms[match(rownames(raw_cts), rownames(tpms)), ]
rownames(tpms) <- rownames(raw_cts)


# Row annot

# add gene symbols
gene_names_df <- data.frame(row.names = rownames(raw_cts))

ens_no_version <- str_remove(rownames(gene_names_df), "\\.\\d+$")
stopifnot(length(ens_no_version) == length(unique(ens_no_version)))

gene_names_df$Symbol <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                              keytype="ENSEMBL", column="SYMBOL", 
                                              multiVals="first")

gene_names_df$Uniq_syms <- scater::uniquifyFeatureNames(rownames(gene_names_df), gene_names_df$Symbol)
gene_names_df$entrez <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                              keytype="ENSEMBL", column="ENTREZID", 
                                              multiVals="first") # there are duplicates in here.

gene_names_df$Gene_name <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                                 keytype="ENSEMBL", column="GENENAME", 
                                                 multiVals="first")

# Sample meta

data_for_DE <- read_tsv(samplesheet) %>%
  as.data.frame() %>%
  dplyr::select(-fq1, -fq2) %>%
  unique() # samplesheet can have more than one row for a given sample (e.g. sequenced on more than one lane)

# samplesheet must have at least sample and group
stopifnot(c("sample", "group") %in% colnames(data_for_DE))

rownames(data_for_DE) <- data_for_DE$sample

# make sure order of samples in the meta data matches the counts
data_for_DE <- data_for_DE[colnames(raw_cts), ]

stopifnot(all(!is.na(data_for_DE$group)))

# Make DDS

stopifnot(identical(rownames(data_for_DE), colnames(raw_cts)))
stopifnot(identical(rownames(gene_names_df), rownames(raw_cts)))
count_data <- list(counts=raw_cts, tpms=tpms[rownames(raw_cts), ])

se <- SummarizedExperiment(assays = count_data, colData = data_for_DE, rowData = gene_names_df)
se <- se[, order(se$group)] # order samples by group

# Add vst
dds <- DESeqDataSet(se, design = ~ group)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

assays(se)$vst <- assay(vsd)

write_rds(se, out_se)

# PCA
sce <- as(se, "SingleCellExperiment")
sce <- runPCA(sce, ntop=5000, ncomponents = 4, exprs_values="vst")

rownames(sce) <- rowData(sce)$Uniq_syms
write_rds(sce, out_sce)


# Output size factors for use with other tools
sizefactors <- 1/dds$sizeFactor
tibble::tibble(sample=names(sizefactors), sizefactor=sizefactors) %>%
  write_tsv(., out_sizeFactors) # "SizeFactors will now contain a factor for every sample which can be used to divide the 4th colun of a bedGraph/bigwig by. Both the aforementioned tools (bamCoverage and genomecov) have options though to directly scale the output by a factor (--scaleFactor or --scale respectively). !! Note though that these options will multiply rather than divide the 4th column with the factor, therefore you would need to provide the reciprocal as mentioned in the code chunk above." https://www.biostars.org/p/413626/#414440

