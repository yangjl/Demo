### Jinliang yang
### 2.17.2015

### input
for(i in 4:6){
  samples = paste0("largedata/sample", i, ".txt")
  
  ### output scripts
  shfile = paste0("largedata/step2_align", i, ".sh")
  slurmfile = paste0("largedata/slurm_step2_align", i, ".sh")
  ### number of CPU to use for mapping
  cpu = 16
  
  
  ######################################################################
  #source("lib/PE_alignment.R")
  setup_PE_alignment(fqfile = samples, shfile = shfile,
                     DBdir="largedata/OS_204_v7/", DBnm="Osative_204_v7", miss=8, cpu=8)
  
  source("lib/setUpslurm.R")
  setUpslurm(slurmsh=slurmfile, wd=NULL, jobid= paste0("align", i),
             codesh= paste("module load gmap/2014-05-15",paste("sh", shfile), sep="\n"))
  ###>>> In this path: cd /home/jolyang/Documents/Github/Demo
  ###>>> [ note: --ntasks=INT, number of cup ]
  ###>>> [ note: --mem=16000, 16G memory ]
  ###>>> RUN: sbatch -p bigmemh --ntasks 16 largedata/slurm_step2_align.sh
  
  
}


###>>> Generate [ largedata/step2_align1.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks 16 largedata/slurm_step2_align1.sh
###>>> Generate [ largedata/step2_align2.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmeml --ntasks 16 largedata/slurm_step2_align2.sh
###>>> Generate [ largedata/step2_align3.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmeml --ntasks 16 largedata/slurm_step2_align3.sh
###>>> Generate [ largedata/step2_align4.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --ntasks 16 largedata/slurm_step2_align4.sh
###>>> Generate [ largedata/step2_align5.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmeml --ntasks 16 largedata/slurm_step2_align5.sh
###>>> Generate [ largedata/step2_align6.sh ] for alignment using GSNAP!
###>>> In this path: cd /home/jolyang/Documents/Github/Demo
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmeml --ntasks 16 largedata/slurm_step2_align6.sh

