library(iSEE)

sce <- readRDS("sce.rds")

rowData(sce) <- rowData(sce)[, setdiff(colnames(rowData(sce)), c("Uniq_syms"))]
search_cols <- rep("", length(colnames(rowData(sce))))
search_cols[grep("\\.rank$", colnames(rowData(sce)))[1]] <- "1 ... 50" # top 50 genes of the first contrast

# To deploy to shinyapp.io, run manually:
# usethis::proj_activate(<<<R_proj_path>>>)
# options(repos = BiocManager::repositories())
# rsconnect::deployApp(appDir=<<<shiny_app_dir>>>, appName = '<<<isee_app_name>>>')
# You may need to follow the instructions at https://docs.posit.co/shinyapps.io/guide/getting_started/#configure-rsconnect to set up your shinyapps.io credentials

iSEE(sce,
     initial=list(ComplexHeatmapPlot(PanelWidth=6L, PanelHeight=1000L,
                                     CustomRows=FALSE,
                                     Assay="vst",
                                     ClusterRows=TRUE,
                                     ShowColumnSelection=FALSE,
                                     ColumnData=c("group"),
                                     AssayCenterRows=TRUE, DataBoxOpen=TRUE,
                                     VisualBoxOpen=TRUE, NamesRowFontSize=14,
                                     NamesColumnFontSize=14, ShowDimNames=c("Rows","Columns"),
                                     RowSelectionSource="RowDataTable1",
                                     LegendPosition="Right",
                                     LegendDirection="Vertical"),
                  RowDataTable(PanelWidth=6L, SearchColumns=search_cols, 
                               HiddenColumns=c("entrez","Gene_name","alias","ens_gene",
                                               grep("\\.LFC$", colnames(rowData(sce)), value = TRUE),
                                               grep("\\.padj$", colnames(rowData(sce)), value = TRUE))),
                  FeatureAssayPlot(DataBoxOpen=TRUE, PanelWidth=6L, Assay="vst",
                                   XAxis="Column data", XAxisColumnData="group",
                                   ColorBy="Column data", ColorByColumnData="group",
                                   FontSize=1.5, PointSize=2),
                  ReducedDimensionPlot(PanelWidth=6L, ColorBy="Column data", ColorByColumnData="group", FontSize=1.5, PointSize=2)
     ),
     appTitle="App for exploring RNA-seq dataset")


