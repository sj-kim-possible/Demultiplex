#!/bin/bash

#SBATCH --partition=bgmp            ### Partition
#SBATCH --job-name=dmux%j           ### Job Name
#SBATCH --output=dmux%j.out         ### File in which to store job output
#SBATCH --error=dmux%j.err          ### File in which to store job error messages
#SBATCH --nodes=1                   ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1           ### Number of cpus per node
#SBATCH --account=bgmp              ### Account used for job submission
#SBATCH --mail-type=ALL             ### Notify user by email when certain events occur
#SBATCH --mail-user=${skim6}@uoregon.edu

read1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
read2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
index1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
index2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
sampleBarcodes="/projects/bgmp/shared/2017_sequencing/indexes.txt"
outDir="/projects/bgmp/skim6/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplexed/"
demultiplexer="/projects/bgmp/skim6/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplex.py"

conda activate bgmp_py310

/usr/bin/time -v $demultiplexer \
    -r1 $read1 \
    -r2 $read2 \
    -i1 $index1 \
    -i2 $index2 \
    -b $sampleBarcodes \
    -o $outDir