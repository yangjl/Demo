### Jinliang Yang
### Edit the ASM465v1.25.gff3


g1 <- read.delim("largedata/OS_indica/Oryza_indica.ASM465v1.25.gff3", header=FALSE, comment="#")

names(g1) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")

g2 <- read.delim("largedata/OS_indica/Osativa_204_v7.0.gene_exons", header=FALSE, comment="#")
names(g2) <- c("seqname", "source", "feature", "start", "end", "score",
               "strand", "frame", "attribute")


g1 <- subset(g1, feature != "repeat_region")
g1 <- g1[!duplicated(g1$attribute),]

g1$feature <- gsub("transcript", "mRNA", g1$feature)

g1$attribute <- gsub("transcript", "mRNA", g1$attribute)
g1$attribute <- gsub("^Name", "ID", g1$attribute)

g1 <- subset(g1, feature %in% c("CDS", "exon", "gene", "mRNA"))


write.table(g1, "largedata/OS_indica/Oryza_indica.ASM465v1.25_edited.gff3", sep="\t",
            row.names=FALSE, quote=FALSE, col.names=FALSE)


g2 <- subset(g1, feature == "mRNA")
tab <- table(g2$attribute)

library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
gr <- import("largedata/OS_indica/Oryza_indica.ASM465v1.25_edited.gff3")
txdb <- makeTranscriptDbFromGFF(file="largedata/OS_indica/Oryza_indica.ASM465v1.25_edited.gff3",
                                format="gff3",
                                dataSource="phytozome",
                                species="oryza")



saveDb(txdb, file="largedata/ASM465v1.25_edited.sqlite")
