### Jinliang yang
### 2.17.2015

cpu <- 8
fq1 <- "SRR1170742.sra_1.fastq"
fq2 <- "SRR1170742.sra_2.fastq"
shfile <- "align1.sh"

fq1 <- "SRR1170744.sra_1.fastq"
fq2 <- "SRR1170744.sra_2.fastq"
shfile <- "align2.sh"

########
source("lib/PE_alignment.R")
setup_PE_alignment(
  shfile = shfile,
  folder = "largedata", cpu = cpu,
  fq1 = fq1,
  fq2 = fq2)

#########
source("lib/setUpslurm.R")
setUpslurm(
  slurmsh = paste0("largedata/run_", shfile),
  codesh = paste0("sh largedata/", shfile),
  oneline=TRUE,
  wd=NULL, jobid= shfile)

