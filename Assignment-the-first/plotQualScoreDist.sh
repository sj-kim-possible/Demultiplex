#!/bin/bash

#SBATCH --partition=bgmp            ### Partition
#SBATCH --job-name=plotDist%j       ### Job Name
#SBATCH --output=plotDist%j.out     ### File in which to store job output
#SBATCH --error=plotDist%j.err      ### File in which to store job error messages
#SBATCH --nodes=1                   ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=4           ### Number of cpus per node
#SBATCH --account=bgmp              ### Account used for job submission
#SBATCH --mail-type=ALL             ### Notify user by email when certain events occur
#SBATCH --mail-user=${skim6}@uoregon.edu

conda activate bgmp_py310

script="/projects/bgmp/skim6/bioinfo/Bi622/Demultiplex/Assignment-the-first/qualScoreDist.py"
read1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
read2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
index1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
index2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"

/usr/bin/time -v $script \
    -f $read1 \
    -l 101 \
    -n 363246735 \
    -p "r1_dist.png"
echo "R1 done"

/usr/bin/time -v $script \
    -f $read2 \
    -l 101 \
    -n 363246735 \
    -p "r2_dist.png"
echo "R2 done"

/usr/bin/time -v $script \
    -f $index1 \
    -l 8 \
    -n 363246735 \
    -p "r3_dist.png"
echo "R3 done"

/usr/bin/time -v $script \
    -f $index2 \
    -l 8 \
    -n 363246735 \
    -p "r4_dist.png"
echo "R4 done"

