#!/usr/bin/env python
#
#.-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#/ / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
#`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
# Quality Score Distribution
#
# Bi622 2022 Demultiplex: Assignment the First
#
# Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
# Average the quality scores at each position for all reads and generate a per nucleotide 
# mean distribution as you did in part 1 of PS4 (in Bi621).
#
# script overview: takes a fastq file and plots the mean distribution per base position
# with a bar plot.
#
#############################################################
#                                                           #
#        Quality Score Distribution per-nucleotide          #
#                                                           #
#############################################################

#### [ port ] ####
import argparse
import numpy as np
import bioinfo
import matplotlib.pyplot as plt
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="A program to plot the distribution of quality scores per base.")
    parser.add_argument("-f", "--file", help="Please specify input *fastq.gz filename", required=True, type = str)
    parser.add_argument("-l", "--readLength", help="Please specify the length of reads", required = True, type=int)
    parser.add_argument("-p", "--plotFilename", help="Please specify the plot filename as .png", required = True, type=str)
    return parser.parse_args()

args = get_args()

#### [ calculate mean distribution ] ####

# init 1D list with 0's
totalQscores = [0]*args.readLength

def populate_list(file: str) -> tuple[list, int]:
    """This function takes an initialized list of base pairs and records the sum of 
    Phred quality scores at each position. Also sums total lines in the file -> this number is records * 4."""
    totalLines = 0
    with gzip.open(file, "rt") as fastq:
        for line in fastq:
            line = line.strip("\n")
            totalLines += 1
            if totalLines % 4 == 0:
                for i in range(len(totalQscores)):
                    totalQscores[i] += bioinfo.convert_phred(line[i])
    return totalQscores, totalLines

# calculate averages
averages, recordSum = populate_list(args.file)
for i in range(len(averages)):
    averages[i] = averages[i]/(recordSum/4)

#############################################################
#                                                           #
#                           plot                            #
#                                                           #
#############################################################

fig, ax = plt.subplots(1, figsize=(15,8))
ax.bar(range(0,args.readLength), averages, color='#15B01A', capsize=3)
plt.xlabel("Read Position")
plt.ylabel("Average Quality Score")
plt.title(f"Per Base Read Quality Distribution using Mean; from file: {args.file}")
plt.savefig(args.plotFilename)
plt.show()
