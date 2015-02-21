### Jinliang Yang
### Feb 20th, 2015

### input
fastqfile = "largedata/sample.txt"

### output scripts
shfile = "largedata/step1_qc.sh"
slurmfile = "largedata/slurm_step1_qc.sh"
### int passes to fastq_quality_filer Minimum quality score to keep
q = 25
### Minimum percent of bases that must have [-q] quality
p = 50

######################################################################
source("lib/PE_qc.R")
PE_qc(fqfile = fastqfile, shfile = shfile, q = q, p = p)

source("lib/setUpslurm.R")
setUpslurm(slurmsh=slurmfile, codesh=paste("sh", shfile), wd=NULL, jobid="qcjob")
