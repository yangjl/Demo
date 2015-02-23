

###################################################
### collect reads
###################################################
ReadCount <- function(eByg = eByg,
                      fastqfile = "largedata/sample.txt"){
  
  
  message(sprintf("###>>> [ %s ] genes loaded ...", length(eByg))) 
  #[1] 39049
  ###############################################################
  
  fq <- read.table(fastqfile, header=TRUE)
  countDF <- data.frame(row.names=names(eByg))
  ldir <- "largedata"
  for(i in 1:nrow(fq)){
    bamfile <- gsub("_1.fastq", ".uniq.bam", fq$fq1[i])
    bamfile <- paste(ldir, bamfile, sep="/")
    aligns <- readGAlignmentsFromBam(bamfile) # Substitute next two lines with this one.
    counts <- countOverlaps(eByg, aligns, ignore.strand=TRUE)
    countDF <- cbind(countDF, counts)
    names(countDF)[ncol(countDF)] <- as.character(fq$fq1[i])
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


