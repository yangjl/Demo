# Jinliang Yang
# updated: 1.17.2015
# setup alignment and collect alignment statistics

########################################################################################
setup_PE_alignment <- function(
  shfile="alignment.sh",
  folder = "largedata", cpu=8,
  fq1 = "SRR1170742.sra_1.fastq",
  fq2 = "SRR1170742.sra_2.fastq"){
  
  cat(paste("# setup alignment", Sys.time(), sep=" "),
      file=shfile, sep="\n")
  
  ###########
  prefix <- gsub("\\..*$", "", fq1)
  uniq <- paste(prefix, "paired_uniq", sep=".")
  bam <- paste(prefix, "paired_uniq.bam", sep=".")
  
  cat(paste("cd", folder),
      #### GSNAP alignment
      # -m miss match
      # -i indel penalty
      # N: look for novel splicing, 1=yes
      # -w: definition of local novel splicing event
      # -A: output sam format
      # -t: number of CUP to use
      # -n: max number of paths to print
      
      paste("gsnap -D ~/dbcenter/OS_indica -d indica_ASM465v1.25_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t", cpu,
            "-n 3 --quality-protocol=sanger --nofails",
            fq1, fq2, "--split-output", prefix, sep=" "),
      ### extract the unique (or reliable) aligned reads
      #http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_FAQ
      paste("samtools view -bS", uniq, ">", bam, sep=" "),
      "",
      file=shfile, sep="\n", append=TRUE)
  system(paste("mv", shfile, folder, sep=" "))
}

### http://research-pub.gene.com/gmap/src/README
#setup_PE_alignment(
#  shfile="align.sh",
#  folder = "largedata", cpu=8,
#  fq1 = "SRR1170742.sra_1.fastq",
#  fq2 = "SRR1170742.sra_2.fastq")


