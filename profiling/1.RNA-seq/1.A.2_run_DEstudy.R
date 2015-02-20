## Jinliang Yang
## 8.11.2014
##

####### quality checking of the reads
library(qrqc)
s.fastq <- readSeqFile(filename="/mnt/02/yangjl/NGS/BD/7-31-14/rep2/Sample_3-2-12-1/3-2-12-1_CAGATC.trimmed.fq",
                       max.length=110)
trimed <- readSeqFile(filename="/mnt/02/yangjl/NGS/BD/4-29-14/Sample_3-2-12-1/3-2-12-1_CAGATC_L006_R1_001.trimmed.fq",
                      max.length=150)
qualPlot(s.fastq)
qualPlot(list("trimmed"=trimed, "untrimmed"=s.fastq))


###################################################
### prepare genomic features
# More Robust: Store Annotations in TranscriptDb
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
library(GenomicFeatures)
txdb <- makeTranscriptDbFromGFF(file="~/dbcenter/OS_indica/Oryza_indica.ASM465v1.25.gff3",
                                format="gff3",
                                dataSource="ftp://ftp.ensemblgenomes.org/pub/plants/release-25/gff3/oryza_indica/Oryza_indica.ASM465v1.25.gff3.gz",
                                species="oryza_indica")
saveDb(txdb, file="cache/Bd192.sqlite")

txdb <- loadDb("cache/Osativa_204_v7.0.sqlite") 
columns(txdb)
keytypes(txdb)

eByg <- exonsBy(txdb, by="gene")
length(eByg)
#[1] 39049

###################################################
### collect read count
###################################################
collect_countDF <- function(bamfile="", eByg=eByg){
  
  countDF <- data.frame(row.names=names(eByg))
  aligns <- readGAlignmentsFromBam(bamfile) # Substitute next two lines with this one.
  counts <- countOverlaps(eByg, aligns, ignore.strand=TRUE)
  countDF <- cbind(countDF, counts)
  names(countDF)[ncol(countDF)] <- bamfile
  return(countDF)  
}

####
countDF1 <- collect_countDF(dir="~/NGS/BD/7-31-14/rep1", eByg=eByg)
countDF2 <- collect_countDF(dir="~/NGS/BD/7-31-14/rep2", eByg=eByg)


repplot <- function(){
  nm <- names(countDF1)
  nm <- gsub("_.*", "", nm)
  par(mfrow=c(3,4))
  for(i in 1:12){
    plot(countDF1[,i], countDF2[,i], xlab="batch 1", ylab="batch 2", main=nm[i])
  }
}
repplot()

#### change the direction of the two WT 
countDF2$tem <- countDF2[, 10]
countDF2[,10] <- countDF2[, 12]
countDF2[,12] <- countDF2[, 15]
countDF2 <- countDF2[, -15]

countDF <- countDF1
for(i in 1:14){
  countDF[,i] <- countDF1[,i] + countDF2[,i]
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

rpkmrep1 <- apply(countDF1, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
rpkmrep2 <- apply(countDF2, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
rpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))


#QC check of the sample reproducibility by computing a correlating matrix and plotting it as a tree.
#Note: the plotMDS function from edgeR is a more robust method for this task.
library(ape)
d <- cor(rpkm, method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)

source("~/Documents/Rcodes/save.append.R")
save.append(list=c("countDF1", "countDF2", "countDF", "rpkmrep1", "rpkmrep2", "rpkm"), 
            file="~/Documents/BDproj/cache/count.RData",
            description="")
