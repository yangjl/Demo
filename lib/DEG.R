## Jinliang Yang
## 8.11.2014
## Using DESeq to analyze differentially expressed genes

deg <- Run_DESeq2(rcdata=as.matrix(rc), colData = targets, designModel = formula(~ rep + tissue))


##########################
DESeq2_pipe <- function(rcdata=rc, colfile="file.txt", 
                        designModel = formula(~ rep + tissue)){
  ################################################
  cols <- read.table(colfile, header=TRUE)
  
  #Raw count data are expected here!
  library(DESeq2)
  ########## WT and WT4C
  dds <- DESeqDataSetFromMatrix(as.matrix(rcdata), colData=cols,
                                design = designModel)
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




