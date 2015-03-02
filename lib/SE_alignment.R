# Jinliang Yang
# updated: 1.17.2015
# setup alignment and collect alignment statistics

setup_SE_alignment <- function(
  fqfile = "largedata/sample.txt", shfile = "largedata/step2_align.sh", 
  DBdir="", DBnm="Osative_204_v7", miss=8, cpu=8
  ){
  ### Note: your fqfile must specify the exact location of your fastq file.
  
  cat(paste("# setup alignment", Sys.time(), sep=" "),
      "#module load gmap/2014-05-15",
      "",
      file=shfile, sep="\n")
  
  fq <- read.table(fqfile, header=TRUE)
  for(i in 1:nrow(fq)){
    ###########
    fq1 <- fq[i, 1] 
    #qc1 <- paste0(fq1, ".qc")
    
    try(prefix <- gsub("\\.fastq$", "", fq1))
    try(prefix <- gsub("\\.fq$", "", fq1))
    uniq <- paste(prefix, "unpaired_uniq", sep=".")
    bam <- paste(prefix, "uniq.bam", sep=".")
    
    cat(#### GSNAP alignment
        # -m miss match
        # -i indel penalty
        # -N: look for novel splicing, 1=yes
        # -w: definition of local novel splicing event
        # -A: output sam format
        # -t: number of CPU to use
        # -n: max number of paths to print
        
        paste("gsnap -D", DBdir,"-d", DBnm, "-m", miss, "-i 2 -N 1 -w 10000 -A sam -t", cpu,
              "-n 3 --quality-protocol=sanger --nofails",
              fq1, "--split-output", prefix, sep=" "),
        ### extract the unique (or reliable) aligned reads
        #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
        paste("samtools view -bS", uniq, ">", bam, sep=" "),
        paste("#rm *concordant* *halfmapping* *nomapping* *paired*"),
        "",
        file=shfile, sep="\n", append=TRUE)
  }
  
  message(sprintf("###>>> Generate [ %s ] for alignment using GSNAP!", shfile))
}


