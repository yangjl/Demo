## Jinliang Yang
## 8.11.2014
## Using DESeq to analyze differentially expressed genes



samples <- "largedata/sample.txt"
design_model <- formula(~ rep + tissue)


################################################
#Raw count data are expected here!
library(DESeq2)

rc <- read.csv("largedata/readcount.csv")
nm <- names(rc)
nm <- gsub("^leaf\\.", "", nm)
nm <- gsub("^root\\.", "", nm)
names(rc) <- nm


targets <- read.table(samples, header=TRUE)
targets$fq1 <- gsub(".*\\/", "", targets$fq1)

##########################
DESeq2_pipe <- function(rcdata=rc, design_table=tgs[c(7:12, 21:26),], 
                        design_model= design_model){
  ########## WT and WT4C
  dds <- DESeqDataSetFromMatrix(as.matrix(rcdata), colData=design_table,
                                design = design_model)
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  res <- as.data.frame(res)
  
  sig <- subset(res, padj < 0.05 & abs(log2FoldChange)>1)
  #myres <- subset(res, padj < 0.05 & abs(log2FoldChange)>2)
  message(sprintf("#>>> DE genes: [ %s ], up: [%s], down: [%s]", nrow(sig), 
                  nrow(subset(sig, log2FoldChange >0)), nrow(subset(sig, log2FoldChange <0)) ))
  return(resOrdered)
}


######## CBF and WT under nonstress
res1 <- DESeq2_pipe(rcdata=as.matrix(rc), design_table= targets, design_model =design_model)
#>>> DE genes: [ 1473 ], up: [840], down: [633]
sig <- subset(res1, padj < 0.05 & abs(log2FoldChange)>1)





