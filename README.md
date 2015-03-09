#RNA-seq on UCD farm

An RNA-seq differential expression analysis pipeline designed for job submission on `UCD farm` cluster using `slurm` system.

# Install

`git clone https://github.com/yangjl/Demo.git`
Remove the largedata folder:
`rmdir largedata`
Then make a link to the `largedata` folder to your directory:
`ln -s /group/jrigrp5/ECL298/Demo/largedata ./largedata`

# Directory

1. `largedata` contains fastq files, sam file, etc.
2. `profiling` R codes
3. `reports` Rmd files to generate reports.

# Tutorial

1. [RNA-seq data](http://rpubs.com/yangjl0930/61344)

2. [RNA-seq Data Analysis](http://rpubs.com/yangjl0930/60157)

# Using pseudo shell to open R
`cd /path/to/Demo`
`srun -p serial --pty --mem 8000 R`
