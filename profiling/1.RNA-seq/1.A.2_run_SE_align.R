### Jinliang yang
### 2.17.2015

### input
source("lib/SE_alignment.R")
source("lib/setUpslurm.R")

######################################################################
#source("lib/PE_alignment.R")
setup_SE_alignment(fqfile = "samples.txt", shfile = "largedata/step2_align.sh",
                   DBdir="largedata/OS_204_v7/", DBnm="Osative_204_v7", miss=8, cpu=32)


setUpslurm(slurmsh="largedata/slurm_step2_align.sh", wd=NULL, jobid= "align",
           codesh= paste("module load gmap/2014-05-15", paste("sh largedata/step2_align.sh"), sep="\n"))

  
