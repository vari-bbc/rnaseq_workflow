args <- commandArgs(trailingOnly = TRUE)

pkges_file <- args[1]

# https://github.com/rstudio/renv/issues/81#issuecomment-497131224
options(repos = BiocManager::repositories()) 

# Get list of packages available on BioConductor
bioc_repos <- BiocManager::repositories()
avail_bioc <- available.packages(repos=bioc_repos[grepl("BioC", names(bioc_repos))])

message("Initiating renv...")
renv::init(bioconductor=TRUE)

# message("Updating packages in renv library...")
# renv::update()

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

            if (pkg %in% avail_bioc[, "Package"]){
                # adding "bioc::" as this has been reported to help prevent installing devel BioC packages and keep packages on the same BioC release. See https://github.com/rstudio/renv/issues/244
                # see also https://rstudio.github.io/renv/reference/install.html for syntax for 'packages' argument
                pkg <- paste0("bioc::", pkg)
            }

            message("Unable to load ", req_pkg, ". Trying to install ", pkg, " ...")
            renv::install(pkg)
		}
}

message("Run renv snapshot...")
renv::snapshot()

sessionInfo()
