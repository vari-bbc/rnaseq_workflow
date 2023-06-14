if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("iSEE", quietly = TRUE))
	BiocManager::install("iSEE")

library(iSEE)

sce <- readRDS("sce.rds")
rowData(sce) <- rowData(sce)[, c("Symbol","entrez","Gene_name")]

iSEE(sce,
     initial=list(FeatureAssayPlot(DataBoxOpen=TRUE, PanelWidth=4L, Assay="vst",
                                   XAxis="Column data", XAxisColumnData="group", 
                                   ColorBy="Column data", ColorByColumnData="group", 
                                   FontSize=1),
                  ReducedDimensionPlot(PanelWidth=3L, ColorBy="Column data", ColorByColumnData="group"),
                  RowDataTable(PanelWidth=5L)
     ), 
     appTitle="App for exploring RNA-seq dataset")


