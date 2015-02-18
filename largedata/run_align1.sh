#!/bin/bash
#SBATCH -D /Users/yangjl/Documents/Github/Demo
#SBATCH -o /Users/yangjl/Documents/Github/Demo/slurm-log/testout-%j.txt
#SBATCH -e /Users/yangjl/Documents/Github/Demo/slurm-log/err-%j.txt
#SBATCH -J align1.sh
set -e
set -u

module load gmap/2014-05-15
sh largedata/align1.sh

python /home/jolyang/bin/send_email.py -s largedata/run_align1.sh
