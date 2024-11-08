args <- commandArgs(trailingOnly = TRUE)

pkges_file <- args[1]

pkges <- read.table(pkges_file)[[1]]

for (pkg in pkges){
        if(!require(pkg, character.only=TRUE)){
            # https://github.com/rstudio/renv/issues/81#issuecomment-497131224
            options(repos = BiocManager::repositories()) 
            renv::install(pkg)
		}
}


renv::snapshot()

