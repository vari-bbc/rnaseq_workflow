args <- commandArgs(trailingOnly = TRUE)

pkges_file <- args[1]

renv::init()
renv::update()

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

renv::settings$bioconductor.version(BiocManager::version())

pkges <- read.table(pkges_file)[[1]]

# https://github.com/rstudio/renv/issues/81#issuecomment-497131224
options(repos = BiocManager::repositories()) 
for (pkg in pkges){
        if(!require(pkg, character.only=TRUE)){
            renv::install(pkg)
		}
}


renv::snapshot()

sessionInfo()
