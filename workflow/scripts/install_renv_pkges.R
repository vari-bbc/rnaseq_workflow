args <- commandArgs(trailingOnly = TRUE)

pkges_file <- args[1]

# https://github.com/rstudio/renv/issues/81#issuecomment-497131224
options(repos = BiocManager::repositories()) 

message("Initiating renv...")
renv::init()

renv::settings$bioconductor.version(BiocManager::version())

message("Updating packages in renv library...")
renv::update()

message("Double check that all dependencies can be loaded.")

pkges <- read.table(pkges_file)[[1]]

for (pkg in pkges){
        req_pkg <- pkg
        if (grepl("\\/", req_pkg)){
            req_pkg <- basename(req_pkg) # remove Github account name if present
            message("Removed Github account name from ", pkg, " to obtain ", req_pkg, ".")
        }
        if (grepl("\\@", req_pkg)){
            req_pkg <- gsub("@\\S+$", "", req_pkg) # remove version number if present
            message("Removed verion number from ", pkg, " to obtain ", req_pkg, ".")
        }
        if(!suppressPackageStartupMessages(require(req_pkg, character.only=TRUE, quietly = TRUE, warn.conflicts = FALSE))){
            message("Unable to load ", req_pkg, ". Trying to install...")
            renv::install(pkg)
		}
}

message("Run renv snapshot...")
renv::snapshot()

sessionInfo()
