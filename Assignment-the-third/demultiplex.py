#!/usr/bin/env python
#
#.-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#/ / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
#`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
# Demultiplex
#
# Bi622 2022 Demultiplex: Assignment the Third
# sj kim
#
# Incorporate feedback from peer code reviews
# Utilize appropriate functions (perhaps you want to import bioinfo???)
# Sufficiently comment your code/use docstrings/use type annotations on functions
# Use unit tests on functions/entire algorithm to ensure it works properly
# Create a useful report for the end user of your code
# Use argparse to "generalize" your code
# Be mindful of "simple" things you can do to optimize your code
# Follow the specifications laid out in Assignment the First for the code
# Unclear? Ask!
#
# Report:
# Percentage of reads from each sample
# Overall amount of index swapping
# Any figures/any other relevant data your code output
#
# script overview: Takes 4 fastq files, hot off the machine, and organizes them into 4 flavors: 
# unknown, low quality, index hopped, and dual matched. Also generates a final report .md 
# to summarize results.

#### [ port ] ####
import argparse
import numpy as np
import bioinfo
import matplotlib.pyplot as plt
import gzip
import math
import itertools

def get_args():
    parser = argparse.ArgumentParser(description="A program to demultiplex raw fastq files.")
    parser.add_argument("-r1", "--bioRead1", help="Please specify biological read 1 input *fastq.gz filename", required=True, type = str)
    parser.add_argument("-r2", "--bioRead2", help="Please specify biological read 2 input *fastq.gz filename. Likely R4.", required = True, type = str)
    parser.add_argument("-i1", "--index1", help="Please specify index 1 input *fastq.gz filename. Likely R2.", required = True, type = str)
    parser.add_argument("-i2", "--index2", help="Please specify index 2 input *fastq.gz filename. Likely R3.", required = True, type = str)
    parser.add_argument("-gz", "--gzip", help="Please specify whether output should be gzipped or not.", action = "store_true")
    parser.add_argument("-b", "--sampleBarcodes", help="Please specify barcodes and sample treatment information.", required = True, type = str)
    parser.add_argument("-o", "--outputDir", help="Please specify output directory for output files.", required = True, type = str)
    return parser.parse_args()

args = get_args()

#############################################################
#                                                           #
#                      Demultiplex                          #
#                                                           #
#############################################################

#### [ high-level functions ] ####
def grab_record(filehandle: str) -> str:
    ''' takes a filehandle and returns a list that holds the record '''
    header = filehandle.readline().strip()
    sequence = filehandle.readline().strip()
    plus = filehandle.readline().strip()
    qscore = filehandle.readline().strip()
    return header, sequence, plus, qscore

def addOneToCount(dictionary: dict, key: str) -> None:
    ''' takes a dictionary and its key and increments one to existing key/value. If not, adds value and sets to 1.'''
    if key in dictionary:
        dictionary[key] += 1
    else:
        dictionary[key] = 1

# make sets of barcodes and reverse complement barcodes for uniqueness
barcodes = set()
revCompBarcodes = set()
with open(args.sampleBarcodes, "r") as indexes:
    indexes.readline() #hacky fix to skip the header- with open is calling readline in the background
    for line in indexes:
        line = line.strip("\n")
        line = line.split("\t")
        barcodes.add(line[-1]) #barcode is the last column of text file
        revCompBarcodes.add(bioinfo.reverse_comp(line[-1]))

flavors = ["unknown", "lowQual", "hopped"]

# writeOut is a dictonary that connects the flavor OR dual-matched indexes to their appropriate file handles
# key: flavor or dual-matched index seq
# value: open filehandle 
writeOut = {}
for flavor in flavors:
    if args.gzip==True:
        r1_name = open(f"{args.outputDir}/{flavor}_R1.fastq.gz", "wt")
        r2_name = open(f"{args.outputDir}/{flavor}_R2.fastq.gz", "wt")
        writeOut[flavor] = [r1_name, r2_name]
    else:
        r1_name = open(f"{args.outputDir}/{flavor}_R1.fastq", "w")
        r2_name = open(f"{args.outputDir}/{flavor}_R2.fastq", "w")
        writeOut[flavor] = [r1_name, r2_name]

for barcode in barcodes:
    if args.gzip==True:
        r1_name = open(f"{args.outputDir}/{barcode}_{barcode}_R1.fastq.gz", "wt")
        r2_name = open(f"{args.outputDir}/{barcode}_{barcode}_R2.fastq.gz", "wt")
        writeOut[barcode] = [r1_name, r2_name]
    else:
        r1_name = open(f"{args.outputDir}/{barcode}_{barcode}_R1.fastq", "w")
        r2_name = open(f"{args.outputDir}/{barcode}_{barcode}_R2.fastq", "w")
        writeOut[barcode] = [r1_name, r2_name]

#### [ stats ] ####
totalRecords = 0
unknownRecords = 0
lowQualityRecords = 0
hoppedRecords = 0
dualMatchedRecords = 0
unkLowCounts = {}
hoppedCounts = {}
dualMatchedCounts = {}


#### [ open all 4 files that need parsing ] ####
with gzip.open(args.bioRead1, "rt") as r1, gzip.open(args.bioRead2, "rt") as r2, gzip.open(args.index1, "rt") as i1, gzip.open(args.index2, "rt") as i2:

    while True:    
        i1_header, i1_seq, i1_plus, i1_qscore = grab_record(i1)
        i2_header, i2_seq, i2_plus, i2_qscore = grab_record(i2)
        r1_header, r1_seq, r1_plus, r1_qscore = grab_record(r1)
        r2_header, r2_seq, r2_plus, r2_qscore = grab_record(r2)

        if i1_header == "":
            break
        
        totalRecords += 1
        
        record1String = f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qscore}\n"
        record2String = f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qscore}\n"

        # index 1 is not a valid barcode
        if i1_seq not in barcodes:
            writeOut["unknown"][0].write(record1String)
            writeOut["unknown"][1].write(record2String)
            unknownRecords += 1
            addOneToCount(unkLowCounts, "unknown")
        
        # index 2 is not a valid barcode
        elif i2_seq not in revCompBarcodes:
            writeOut["unknown"][0].write(record1String)
            writeOut["unknown"][1].write(record2String)
            unknownRecords += 1
            addOneToCount(unkLowCounts, "unknown")

        # if the average quality score of index 1 is below 30
        elif bioinfo.qual_score(i1_qscore) < 30:
            writeOut["lowQual"][0].write(record1String)
            writeOut["lowQual"][1].write(record2String)
            lowQualityRecords += 1
            addOneToCount(unkLowCounts, "lowQual")
        
        # if avg qual score of index 2 is below 30
        elif bioinfo.qual_score(i2_qscore) < 30:
            writeOut["lowQual"][0].write(record1String)
            writeOut["lowQual"][1].write(record2String)
            lowQualityRecords += 1
            addOneToCount(unkLowCounts, "lowQual")
        
        # check hopped
        elif i1_seq != bioinfo.reverse_comp(i2_seq):
            writeOut["hopped"][0].write(record1String)
            writeOut["hopped"][1].write(record2String)
            hoppedRecords += 1
            addOneToCount(hoppedCounts, f"{i1_seq}-{bioinfo.reverse_comp(i2_seq)}")

        # dual matched
        elif i1_seq == bioinfo.reverse_comp(i2_seq):
            r1_header += f"_{i1_seq}_{bioinfo.reverse_comp(i2_seq)}"
            writeOut[i1_seq][0].write(f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qscore}\n")
            r2_header += f"_{i1_seq}_{bioinfo.reverse_comp(i2_seq)}"
            writeOut[i1_seq][1].write(f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qscore}\n")
            dualMatchedRecords += 1
            addOneToCount(dualMatchedCounts, f"{i1_seq}-{bioinfo.reverse_comp(i2_seq)}")
        
# close the write out files
for filename in writeOut:
    writeOut[filename][0].close()
    writeOut[filename][1].close()

# compile stats into .txt
with open(f"{args.outputDir}/summary_stats.txt", "w") as stats:
    stats.write(f"Total Number of Records: {totalRecords}\n")
    stats.write(f"Unknown Records: {unknownRecords}, {unknownRecords/totalRecords*100}%\n")
    stats.write(f"Low Quality Records: {lowQualityRecords}, {lowQualityRecords/totalRecords*100}%\n")
    stats.write(f"Index Hopped Records: {hoppedRecords}, {hoppedRecords/totalRecords*100}%\n")
    stats.write(f"Dual-Matched Records: {dualMatchedRecords}, {dualMatchedRecords/totalRecords*100}%\n")
    stats.write(f"Unknown and Low Quality Totals:\n")
    stats.write(f"Index1-Index2\tNumber of Records\tPercent of Total Records\n")
    for record in unkLowCounts:
        stats.write(f"{record}\t{unkLowCounts[record]}\t{unkLowCounts[record]/totalRecords}\n")
    stats.write(f"Hopped Totals:\n")
    stats.write(f"Index1-Index2\tNumber of Records\tPercent of Total Records\n")
    for record in hoppedCounts:
        stats.write(f"{record}\t{hoppedCounts[record]}\t{hoppedCounts[record]/totalRecords}\n")
    stats.write(f"Dual-Matched Totals:\n")
    stats.write(f"Index1-Index2\tNumber of Records\tPercent of Total Records\n")
    for record in dualMatchedCounts:
        stats.write(f"{record}\t{dualMatchedCounts[record]}\t{dualMatchedCounts[record]/totalRecords}\n")
    