## Jinliang Yang
## 8.11.2014
##


featureDB <- "largedata/Osativa_204_v7.0.sqlite"
samples <- "largedata/sample.txt"


#############################
#samtools view -h -o out.sam in.bam
#samtools view -h leaf.rep1_1.fastq.concordant_uniq.bam -o leaf.rep1_1.fastq.concordant_uniq.sam
#htseq-count leaf/leaf.rep3.concordant_uniq OS_indica/Oryza_indica.ASM465v1.25.gff3 -i Parent > test2.out


source("lib/readcount.R")
library(GenomicFeatures)
library(GenomicAlignments)
#### load db
txdb <- loadDb(featureDB) 
#columns(txdb)
#keytypes(txdb)
eByg <- exonsBy(txdb, by="gene")

rc <- ReadCount(eByg = eByg, fastqfile = samples)

write.table(rc, "largedata/readcount.csv", sep=",", row.names=FALSE, quote=FALSE)


rc <- read.csv("largedata/readcount.csv")
rpkm <- apply(rc, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
write.table(rpkm, "largedata/rpkm.csv", sep=",", row.names=FALSE, quote=FALSE)


#rpkm <- read.csv("largedata/rpkm.csv")
#tem <- rpkm[,1]
#plot(density(tem[tem>10]))

