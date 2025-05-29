args <- commandArgs(trailingOnly = TRUE)

gtf_file <- args[1]
orgdb <- args[2]
renv_rproj_dir <- args[3]
out_se <- args[4]
out_sce <- args[5]
out_sizeFactors <- args[6]
out_txi <- args[7]
out_strand <- args[8]

# use the packages in the renv
renv::load(renv_rproj_dir)

star_dir <- "results/star"
salmon_dir <- "results/salmon"
samplesheet <- "config/samplesheet/units.tsv"

# Packages loaded

library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(edgeR)
library(rtracklayer)
# load the org.db for your organism
library(orgdb, character.only=TRUE)
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

message(str_glue("Salmon-inferred library type is {lib_type}."))

fileConn <- file(out_strand)
writeLines(lib_type, con = fileConn)
close(fileConn)

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

# extract the first transcript from one of the salmon result files
test_gtf_tx <- stringr::str_split_1(readLines(files[1], 2)[2], "\\t")[1]
message("Using ", test_gtf_tx, " to decide whether to turn on ignoreTxVersion in tximport.")

tximport_ignore_tx_ver <- NULL
tximport_ignore_tx_ver <- if (test_gtf_tx %in% tx2gene$TXNAME) {
  tximport_ignore_tx_ver <- FALSE
} else if (str_remove(test_gtf_tx, "\\.\\d+$") %in% tx2gene$TXNAME) {
  tximport_ignore_tx_ver <- TRUE
} else{
  stop("No match between GTF and Salmon transcript names with or without removing transcript version.")
}
message("The 'ignoreTxVersion' parameter in tximport will be set to ", as.character(tximport_ignore_tx_ver))

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion=tximport_ignore_tx_ver)
write_rds(txi.salmon, out_txi)

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

## add alias column to gene_names_df
gtf_df <- rtracklayer::import(gtf_file) |> 
  as.data.frame() |>
  dplyr::select(gene_id, gene_name) |>
  dplyr::distinct() |>
  column_to_rownames("gene_id")

stopifnot(all(rownames(gtf_df) %in% rownames((gene_names_df))))
gene_names_df$alias <- gtf_df$gene_name[match(rownames(gene_names_df), rownames(gtf_df))]

# Sample meta

data_for_DE <- read_tsv(samplesheet) %>%
  as.data.frame() %>%
  dplyr::select(-fq1, -fq2) %>%
  unique() # samplesheet can have more than one row for a given sample (e.g. sequenced on more than one lane)

stopifnot(length(data_for_DE$sample) == length(unique(data_for_DE$sample)))

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

rowData(sce)$ens_gene <- rownames(sce)
rownames(sce) <- rowData(sce)$Uniq_syms
write_rds(sce, out_sce)


# Output size factors for use with other tools
sizefactors <- 1/dds$sizeFactor
tibble::tibble(sample=names(sizefactors), sizefactor=sizefactors) %>%
  write_tsv(., out_sizeFactors) # "SizeFactors will now contain a factor for every sample which can be used to divide the 4th colun of a bedGraph/bigwig by. Both the aforementioned tools (bamCoverage and genomecov) have options though to directly scale the output by a factor (--scaleFactor or --scale respectively). !! Note though that these options will multiply rather than divide the 4th column with the factor, therefore you would need to provide the reciprocal as mentioned in the code chunk above." https://www.biostars.org/p/413626/#414440

sessionInfo()
