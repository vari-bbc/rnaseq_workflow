options(usethis.allow_nested_project = TRUE)

args <- commandArgs(trailingOnly = TRUE)

R_proj_name <- args[1]

for (pkg in c("usethis","renv","BiocManager")){
	if(!require(pkg, character.only=TRUE)){
        # specify 'repos' to avoid being asked to choose mirror.
    	install.packages(pkg, repos = "https://cloud.r-project.org/")
	}
}

usethis::create_project(path = R_proj_name, open = TRUE, rstudio = TRUE)

sessionInfo()

