library(optparse, quietly = T)
library(rtracklayer, quietly = T)
library(tidyverse, quietly = T)
library(data.table, quietly = T)
library(Biostrings, quietly = T)

option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
                help="path to raw gff and peptide data (downloaded from ***.sh)"),
  make_option(c("-o", "--output"), default=NULL,
                help="path to output for parsed data. default is input_dir/annotated/"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## function to get filepaths for raw data
## modified from GENESPACE
## https://github.com/jtlovell/GENESPACE/blob/390341499ee1d2ccd5e1a894c4bd7c1bd20a3dda/R/init_genespace.R#L231-L273
get_filePaths <- function(path, pattern){
  if(!dir.exists(path))
    stop("can't find the directory!")
  
  pattern <- as.character(pattern)
  ps <- sapply(path, USE.NAMES=T, simplify=F, function(x)
    list.files(
      path=x,
      pattern=pattern,
      full.names=T,
      recursive=T))
  ps0 <- unlist(ps)
  return(ps0)
}

## function to parse NCBI genome into annotation
## modified from GENESPACE
## https://github.com/jtlovell/GENESPACE/blob/master/R/parse_annotations.R
parse_annot <- function(gffIn, gffOut, pepIn, pepOut){
  
  # read in gene gff
  gff <- data.table(data.frame(rtracklayer::readGFF(
    filepath=gffIn,
    filter=list(type="gene"),
    tags=c("gene","gene_biotype"))))
  gff <- subset(gff, gene_biotype=="protein_coding")
  gff <- subset(gff, !duplicated(gene))
  gff[,`:=`(seqid = as.character(seqid), gene = as.character(gene))]
  gchr <- gff[,list(nGene = .N), by = "seqid"]
  
  # read in region gff
  chrIDs <- data.table(data.frame(rtracklayer::readGFF(
    filepath = gffIn,
    filter = list(type = "region"),
    tags = c("chromosome"))))
  chrIDs[,`:=`(seqid = as.character(seqid), chromosome = as.character(chromosome))]
  
  uchrs <- c(unique(chrIDs$chromosome), unique(chrIDs$seqid))
  uchrs <- uchrs[!duplicated(uchrs)]
  
  # choose largest regions / chromosome as rep
  gid <- chrIDs[,list(nbp = end - start), by = c("chromosome","seqid")]
  gchr <- merge(gchr, gid, by = "seqid", all.x = T)
  gchr[,isBest := nGene == max(nGene) & !is.na(chromosome) & chromosome != "Unknown", by = "chromosome"]
  gchr[,chr := ifelse(isBest, chromosome, seqid)]
  setorder(gchr, -nbp)
  gchr <- subset(gchr, !duplicated(seqid))
  
  # rename sequences in gene gff with regions
  gff <- merge(gff, gchr[,c("seqid", "chr")], by = "seqid")
  uchrs <- uchrs[uchrs %in% unique(gff$chr)]
  
  # order genes
  gff[,chr := factor(chr, levels = uchrs)]
  gff$chr[is.na(gff$chr)] <- gff$seqid[is.na(gff$chr)]
  setorder(gff, chr, start, end, na.last = T)
  
  # read in peptide fa
  fa <- readAAStringSet(pepIn)
  
  # rename peptide headers
  tmp <- sapply(names(fa), function(x) strsplit(x, " ")[[1]][2])
  tmp2 <- substr(tmp, 7, nchar(tmp)-1)
  names(fa) <- tmp2
  nfa <- length(fa)
  
  # drop duplicates, keeping the longest
  fa <- fa[order(-width(fa))]
  fa <- fa[!duplicated(names(fa))]
  fa <- AAStringSet(gsub(".","",fa, fixed = T))
  fa <- fa[width(fa) > 20]
  
  # join the two
  gff <- subset(gff, gene %in% names(fa))
  fa <- fa[gff$gene]
  
  if(nrow(gff) == 0){
    stop("PARSING FAILED - is this a ncbi-formatted annotation?\n")
  }else{
    
    # write output
    gff <- with(gff, data.table(
      chr = chr,
      start = start,
      end = end,
      id = gene,
      strand = strand,
      ord = 1:length(gene)))
    
    fwrite(gff, file = gffOut, sep = "\t", quote = F)
    writeXStringSet(fa, filepath = pepOut)
  }
}


## find species to parse
species <- list.files(opt$directory)

## make directory for output
output <- opt$output
if (is.null(output)){output <- file.path(opt$directory, 'annotated')}
if (!dir.exists(output)){dir.create(output)}



## ID inputs and outputs
gffIn <- get_filePaths(opt$directory, "gff")
pepIn <- get_filePaths(opt$directory, "faa")

for (i in 1:length(gffIn)){
  gIN <- gffIn[i] # path to gene
  pIN <- pepIn[i] # path to protein
  
  gB <- basename(gIN)#, opt$directory) # file ID
  pB <- basename(pIN)#, opt$directory) # file ID
  
  gOUT <- paste(output, gB, sep="/") # output path
  pOUT <- paste(output, pB, sep="/") # output path
  
  gIN <- as.character(gIN)
  gOUT <- as.character(gOUT)
  pIN <- as.character(pIN)
  pOUT <- as.character(pOUT)
  parse_annot(gIN, gOUT, pIN, pOUT)
}
