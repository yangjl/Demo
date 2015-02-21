# Jinliang Yang
# updated: 1.17.2015
# setup alignment and collect alignment statistics

setup_PE_alignment <- function(
  fqfile = "largedata/sample.txt", shfile = "largedata/step2_align.sh", cpu=8
  ){
  
  cat(paste("# setup alignment", Sys.time(), sep=" "),
      "#module load gmap/2014-05-15",
      "",
      file=shfile, sep="\n")
  
  fq <- read.table(fqfile, header=TRUE)
  for(i in 1:nrow(fq)){
    ###########
    fq1 <- paste0("largedata/", fq$fq1[i]) 
    #qc1 <- paste0(fq1, ".qc")
    fq2 <- paste0("largedata/", fq$fq2[i]) 
    #qc2 <- paste0(fq2, ".qc")
    
    prefix <- gsub("_1\\.fastq$", "", fq1)
    uniq <- paste(prefix, "concordant_uniq", sep=".")
    bam <- paste(prefix, "concordant_uniq.bam", sep=".")
    
    cat(#### GSNAP alignment
        # -m miss match
        # -i indel penalty
        # -N: look for novel splicing, 1=yes
        # -w: definition of local novel splicing event
        # -A: output sam format
        # -t: number of CPU to use
        # -n: max number of paths to print
        
        paste("gsnap -D largedata/OS_indica -d ASM465v1.25_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t", cpu,
              "-n 3 --quality-protocol=sanger --nofails",
              fq1, fq2, "--split-output", prefix, sep=" "),
        ### extract the unique (or reliable) aligned reads
        #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
        paste("samtools view -bS", uniq, ">", bam, sep=" "),
        "",
        file=shfile, sep="\n", append=TRUE)
  }
  
  message(sprintf("###>>> Generate [ %s ] for alignment using GSNAP!", shfile))
}


