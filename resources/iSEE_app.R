library(iSEE)

sce <- readRDS("sce.rds")

# sort by significance of the first occurring contrast
sce <- sce[order(rowData(sce)[[grep("padj$", colnames(rowData(sce)))[1]]]), ]
rowData(sce)$row_number <- 1:nrow(sce)

rowData(sce) <- rowData(sce)[, setdiff(colnames(rowData(sce)), c("Uniq_syms"))]
search_cols <- rep("", length(colnames(rowData(sce))))
search_cols[which(colnames(rowData(sce)) == "row_number")] <- "0 ... 50"

# To deploy to shinyapp.io, run manually:
# usethis::proj_activate(<<<R_proj_path>>>)
# options(repos = BiocManager::repositories())
# rsconnect::deployApp(appDir=<<<shiny_app_dir>>>, appName = '<<<isee_app_name>>>')
# You may need to follow the instructions at https://docs.posit.co/shinyapps.io/guide/getting_started/#configure-rsconnect to set up your shinyapps.io credentials

iSEE(sce,
     initial=list(ComplexHeatmapPlot(PanelWidth=5L, PanelHeight=1000L,
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
                  RowDataTable(PanelWidth=5L, SearchColumns=search_cols, HiddenColumns=c("entrez","Gene_name","alias","ens_gene","row_number")),
                  FeatureAssayPlot(DataBoxOpen=TRUE, PanelWidth=4L, Assay="vst",
                                   XAxis="Column data", XAxisColumnData="group",
                                   ColorBy="Column data", ColorByColumnData="group",
                                   FontSize=1.5, PointSize=2),
                  ReducedDimensionPlot(PanelWidth=3L, ColorBy="Column data", ColorByColumnData="group", FontSize=1.5, PointSize=2)
     ),
     appTitle="App for exploring RNA-seq dataset")


