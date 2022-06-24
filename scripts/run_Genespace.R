# load libraries
library(GENESPACE)
library(optparse)

option_list = list(
  make_option(c("-g", "--genomes"), type="character", default=NULL,
              help="path to raw genome directory"),
  make_option(c("-s", "--species"), type="character", default=NULL,
              help="comma-separated list of species IDs (first will be scaffold for pan-genome)", metavar="character"),
  make_option(c("-w", "--workingdir"), type="character", default=NULL,
              help="path to working directory [default= %default]", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (is.null(opt$workingdir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# path to working directory (will create)
runwd<-file.path(opt$workingdir)

# initialize run
gids<-unlist(strsplit(opt$species, split=","))

gpar<-init_genespace(
  genomeIDs=gids,     # vector of unique names to call genomes
  speciesIDs=gids,    # subdirectories of rawGenomeDir containing data for each genomeID
  versionIDs=gids,    # subdirectories of rawGenomeDir/speciesIDs pointing to version of data
  ploidy=rep(1,length(gids)),    # integer string specifying ploidy of genome assemblies(?)
  wd=runwd,           # where output is stored
  gffString="gff",    # string to use to search for gff files
  pepString="pep",    # string to use to search for pep files
  path2orthofinder="orthofinder",
  path2diamond="diamond",
  path2mcscanx="scripts/MCScanX",
  rawGenomeDir=opt$genomes)

# format raw annotations
parse_ncbi(
  gsParam=gpar)

# run orthofinder
gpar<-run_orthofinder(gsParam=gpar, overwrite=FALSE)

# synteny search
gpar<-synteny(gsParam=gpar)

# pangenome annotation
pg<-pangenome(gpar)
