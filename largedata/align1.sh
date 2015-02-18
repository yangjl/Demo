# setup alignment 2015-02-17 23:13:40
cd largedata
gsnap -D ~/dbcenter/OS_indica -d indica_ASM465v1.25_gsnap -m 10 -i 2 -N 1 -w 10000 -A sam -t 8 -n 3 --quality-protocol=sanger --nofails SRR1170742.sra_1.fastq SRR1170742.sra_2.fastq --split-output SRR1170742
samtools view -bS SRR1170742.paired_uniq > SRR1170742.paired_uniq.bam

