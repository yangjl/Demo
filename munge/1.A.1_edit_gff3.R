### Jinliang Yang
### Edit the ASM465v1.25.gff3


g1 <- read.delim("largedata/OS_indica/Oryza_indica.ASM465v1.25.gff3", header=FALSE, comment="#")

names(g1) <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute")

g2 <- read.delim("largedata/OS_indica/Osativa_204_v7.0.gene_exons", header=FALSE, comment="#")
names(g2) <- c("seqname", "source", "feature", "start", "end", "score",
               "strand", "frame", "attribute")


g1 <- subset(g1, feature != "repeat_region")
g1$feature <- gsub("transcript", "mRNA", g1$feature)

g1$attribute <- gsub("transcript", "mRNA", g1$attribute) 
write.table(g1, "largedata/OS_indica/Oryza_indica.ASM465v1.25_edited.gff3", sep="\t",
            row.names=FALSE, quote=FALSE, col.names=FALSE)