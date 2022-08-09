#!/bin/bash

#SBATCH --partition=bgmp            ### Partition
#SBATCH --job-name=pigz%j           ### Job Name
#SBATCH --output=pigz%j.out         ### File in which to store job output
#SBATCH --error=pigz%j.err          ### File in which to store job error messages
#SBATCH --nodes=1                   ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=12          ### Number of cpus per node
#SBATCH --account=bgmp              ### Account used for job submission
#SBATCH --mail-type=ALL             ### Notify user by email when certain events occur
#SBATCH --mail-user=${skim6}@uoregon.edu

unzipped="/projects/bgmp/skim6/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplexed/"

/usr/bin/time -v pigz $unzipped/*.fastq