#RNA-seq Experiment

# Install

`git clone https://github.com/yangjl/Demo.git`
In addition, you have to manually copy the `largedata` folder to your directory:
`cp -r /group/grigrp5/ECL298/largedata/ .`


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
