#!/usr/bin/env bash

#SBATCH --job-name="Demultiplex"
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --time=5-0

conda activate bgmp_py37

/usr/bin/time \
/Users/agray11/bioinformatics/Bi622/Demultiplexing/demultiplexing-alyxgray7/Assignment-the-third/demux.py \
-r1 /Users/agray11/bioinformatics/Bi622/Demultiplexing/demultiplexing-alyxgray7/Assignment-the-first/test_R1.fastq \
-r2 /Users/agray11/bioinformatics/Bi622/Demultiplexing/demultiplexing-alyxgray7/Assignment-the-first/test_R2.fastq \
-r3 /Users/agray11/bioinformatics/Bi622/Demultiplexing/demultiplexing-alyxgray7/Assignment-the-first/test_R3.fastq \
-r4 /Users/agray11/bioinformatics/Bi622/Demultiplexing/demultiplexing-alyxgray7/Assignment-the-first/test_R4.fastq \
-cutoff 20

/usr/bin/time \
gzip /projects/bgmp/agray11/bioinfo/Bi622/demultiplexing-alyxgray7/Assignment-the-third/output/*.fastq