## Jinliang Yang
## 8.11.2014
##

samtools view -h -o out.sam in.bam
samtools view -h leaf.rep1_1.fastq.concordant_uniq.bam -o leaf.rep1_1.fastq.concordant_uniq.sam
htseq-count leaf/leaf.rep3.concordant_uniq OS_indica/Oryza_indica.ASM465v1.25.gff3 -i Parent > test2.out




library(GenomicFeatures)
library(GenomicAlignments)

txdb <- makeTranscriptDbFromGFF(file="largedata/OS_indica/Oryza_indica.ASM465v1.25_edited.gff3",
                                format="gff3",
                                dataSource="ftp://ftp.ensemblgenomes.org/pub/plants/release-25/gff3/oryza_indica/Oryza_indica.ASM465v1.25.gff3.gz",
                                species="oryza_indica")

saveDb(txdb, file="cache/Bd192.sqlite")

txdb <- loadDb("largedata/Osativa_204_v7.0.sqlite") 
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
countDF1 <- collect_countDF(bamfile="largedata/SRR1170742.concordant_uniq.bam", eByg=eByg)
countDF2 <- collect_countDF(bamfile="largedata/SRR1170744.concordant_uniq.bam", eByg=eByg)


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
