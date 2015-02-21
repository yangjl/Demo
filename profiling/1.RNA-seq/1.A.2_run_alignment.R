### Jinliang yang
### 2.17.2015

### input
fastqfile = "largedata/sample.txt"

### output scripts
shfile = "largedata/step2_align.sh"
slurmfile = "largedata/slurm_step2_align.sh"
### number of CPU to use for mapping
cpu = 8


######################################################################
source("lib/PE_alignment.R")
setup_PE_alignment(fqfile = fastqfile, shfile = shfile, cpu=cpu)

source("lib/setUpslurm.R")
setUpslurm(slurmsh=slurmfile, wd=NULL, jobid="mapping",
           codesh= paste("module load gmap/2014-05-15",paste("sh", shfile), sep="\n"))

gsnap -D ~/Documents/Github/Demo/largedata/OS_indica -d ASM465v1.25_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t 8 -n 3 --quality-protocol=sanger --nofails largedata/leaf/leaf.rep1_1.fastq largedata/leaf/leaf.rep1_2.fastq --split-output largedata/leaf/leaf.rep1_1.fastq