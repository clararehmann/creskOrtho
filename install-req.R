# make local library (if using Talapas)
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

rpacks <- c('dplyr', 'patchwork', 'optparse', 'Seurat', 'tidyverse', 'data.table')
rinstall <- function(package){
    install.packages(package, repos='http://cran.us.r-project.org')
}

for (pack in rpacks){
    tryCatch(find.package(pack), error=function(e) rinstall(pack))
}

bpacks <- c('Biostrings', 'rtracklayer')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
binstall <- function(package){
    BiocManager::install(package)
}

for (pack in bpacks){
    tryCatch(find.package(pack), error=function(e) binstall(pack))
}
