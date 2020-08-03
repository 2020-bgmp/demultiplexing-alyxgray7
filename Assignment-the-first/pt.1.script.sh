#!/usr/bin/env bash

#SBATCH --job-name="R4.qscore_distributions"
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=7
#SBATCH --time=5-0

conda activate bgmp_py37
/usr/bin/time \
/projects/bgmp/agray11/bioinfo/Bi622/demultiplexing-alyxgray7/Assignment-the-first/distributions.py \
-f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-l 101 \
-o R4 

