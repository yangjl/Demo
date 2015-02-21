# Jinliang Yang
# updated: 1.17.2015

PE_qc <- function(fqfile = "largedata/sample.txt", shfile = "largedata/step1_qc.sh",
                  q = 25, p =50){
  
  fq <- read.table(fqfile, header=TRUE)
  
  cat(paste("#quality checking", Sys.time(), sep=" "),
      paste("module load FASTX-Toolkit/0.0.13.2-goolf-1.4.10"),
      "",
      file=shfile, sep="\n")
  
  for(i in 1:nrow(fq)){
    fq1 <- paste0("largedata/", fq$fq1[i])
    stat1 <- paste(fq1, "stat", sep=".")
    qcpng1 <- paste(fq1, "qc.png", sep=".")
    nudpng1 <- paste(fq1, "nud.png", sep=".")
    gz1 <- paste(fq1, "qc", sep=".")
    res1 <- paste(fq1, "res", sep=".")
    
    fq2 <- paste0("largedata/", fq$fq2[i])
    stat2 <- paste(fq2, "stat", sep=".")
    qcpng2 <- paste(fq2, "qc.png", sep=".")
    nudpng2 <- paste(fq2, "nud.png", sep=".")
    gz2 <- paste(fq2, "qc", sep=".")
    res2 <- paste(fq2, "res", sep=".")
    
    cat(
      #### quality checking stat
      #fastx_quality_stats -Q33 -i leaf.rep3_1.fastq -o leaf.rep3_1.stat2 
      paste("fastx_quality_stats -Q33 -i", fq1, "-o", stat1, sep=" "),
      paste("fastx_quality_stats -Q33 -i", fq2, "-o", stat2, sep=" "),
      
      #fastq_quality_boxplot_graph.sh -i leaf.rep3_1.stat2 -o leaf.rep3_1.qc.png -t "leaf.rep3_1"
      paste("fastq_quality_boxplot_graph.sh -i", stat1, "-o", qcpng1, "-t", fq1),
      paste("fastq_quality_boxplot_graph.sh -i", stat2, "-o", qcpng2, "-t", fq2),
      #fastx_nucleotide_distribution_graph.sh -i leaf.rep3_1.stat1 -o leaf.rep3_1.nud.png -t "leaf.rep3_1"
      paste("fastx_nucleotide_distribution_graph.sh -i", stat1, "-o", nudpng1, "-t", fq1),
      paste("fastx_nucleotide_distribution_graph.sh -i", stat2, "-o", nudpng2, "-t", fq2),
      
      #### quality filter
      ### Usually qualities over 30 are considered good and under 20 are considered bad
      #fastq_quality_filter -Q33 -q 25 -p 50 -i leaf.rep3_1.fastq -o leaf.rep3_1_qf.fq -v > leaf.rep3_1.qf
      paste("fastq_quality_filter -Q33 -q", q, "-p", p, "-z -i", fq1, "-o", gz1, "-v >", res1),
      paste("fastq_quality_filter -Q33 -q", q, "-p", p, "-z -i", fq2, "-o", gz2, "-v >", res2),
      "",
      file=shfile, sep="\n", append=TRUE)
    
  }
  message(sprintf("###>>> Generate [ %s ] for quality checking using FASTX!", shfile))
}



