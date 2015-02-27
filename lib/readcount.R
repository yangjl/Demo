
###################################################
### collect reads
###################################################
ReadCount <- function(featureDB = "largedata/ASM465v1.25_edited.sqlite",
                      bamfile = "largedata/sample.txt"){
  library(GenomicFeatures)
  library(GenomicAlignments)
  
  #### load db
  txdb <- loadDb(featureDB) 
  #columns(txdb)
  #keytypes(txdb)
  eByg <- exonsBy(txdb, by="gene")
  
  message(sprintf("###>>> loaded [ %s ] genes for reads counting <<<###", length(eByg))) 
  #[1] 39049
  ###############################################################
  
  fq <- read.table(bamfile, header=FALSE)
  countDF <- data.frame(row.names=names(eByg))
  #ldir <- "largedata"
  for(i in 1:nrow(fq)){
    bam <- as.character(fq$V1[i])
    message(sprintf("###>>> start counting for [ %s ]", bam))
    
    #bamfile <- paste(ldir, bamfile, sep="/")
    aligns <- readGAlignmentsFromBam(bam) # Substitute next two lines with this one.
    counts <- countOverlaps(eByg, aligns, ignore.strand=TRUE)
    countDF <- cbind(countDF, counts)
    names(countDF)[ncol(countDF)] <- bam
  }
  message("### DONE!")
  return(countDF)
}


###################################################
### Compute RPKM
###################################################
returnRPKM <- function(counts, gffsub=eByg) {
  geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
  millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
  rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
  rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
  return(rpkm)
}


