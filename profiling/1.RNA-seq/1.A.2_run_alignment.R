### Jinliang yang
### 2.17.2015

### input
samples = "largedata/sample.txt"

### output scripts
shfile = "largedata/step2_align.sh"
slurmfile = "largedata/slurm_step2_align.sh"
### number of CPU to use for mapping
cpu = 16


######################################################################
source("lib/PE_alignment.R")
setup_PE_alignment(fqfile = samples, shfile = shfile, cpu=cpu)

source("lib/setUpslurm.R")
setUpslurm(slurmsh=slurmfile, wd=NULL, jobid="mapping",
           codesh= paste("module load gmap/2014-05-15",paste("sh", shfile), sep="\n"))
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks 16 largedata/slurm_step2_align.sh
