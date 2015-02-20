# Jinliang Yang
# function to prepare slurm script

setUpslurm <- function(slurmsh="largedata/GenSel/CL_test.sh",
                       oneline=TRUE,
                       codesh="myscript.sh",
                       wd=NULL, jobid="myjob"
                       ){
    
    #message(sprintf("###>>> cp from Introgression, tailored for pvpDiallel"))  
  
    ##### setup working directory
    if(is.null(wd)){
       wd <- getwd()
    }
    dir.create("slurm-log", showWarnings = FALSE)
    sbath <- paste0(wd, "/slurm-log/")
    sbatho <- paste0(sbath, "testout-%j.txt")
    sbathe <- paste0(sbath, "err-%j.txt")
    sbathJ <- jobid
    
    #### parameters pass to slurm script
    cat(paste("#!/bin/bash"),
        #-D sets your project directory.
        #-o sets where standard output (of your batch script) goes.
        #-e sets where standard error (of your batch script) goes.
        #-J sets the job name.
        paste("#SBATCH -D", wd, sep=" "),
        paste("#SBATCH -o", sbatho, sep=" "),
        paste("#SBATCH -e", sbathe, sep=" "),
        paste("#SBATCH -J", sbathJ, sep=" "),
        "set -e",
        "set -u",
        "",
        "module load gmap/2014-05-15",
        file=slurmsh, sep="\n", append=FALSE);
    
    #### attach some sh scripts
    if(oneline){
       cat(codesh, file=slurmsh, sep="\n", append=TRUE) 
    }else{
        cat(paste("sh", codesh),
            file=slurmsh, sep="\n", append=TRUE)
    }
    
    #### warning and message
    cat("",
        paste("#python /home/jolyang/bin/send_email.py -s", slurmsh),
        file=slurmsh, sep="\n", append=TRUE);
    
    message(paste("###>>> In this path: cd ", wd, sep=""), "\n",
            paste("###>>> [ note: --ntasks=INT, number of cup ]"),"\n",
            paste("###>>> RUN: sbatch -p bigmemh --ntasks 8", slurmsh),
            "")
    
}

