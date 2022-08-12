#!/usr/bin/env python

#script overview: takes the index hopped output from summary_stats.txt (hopped_output.txt) and sorts the 
# hopped indexes by order of frequency from most to least in a new text file.

import matplotlib.pyplot as plt
import numpy as np

hopped = {}
sortedHopped = {}

with open("hopped_output.txt", "r") as hopped_infile:
    hopped_infile.readline() #skip header
    for line in hopped_infile:
        line = line.strip()
        line = line.split("\t")
        hopped[line[0]] = (int(line[1]), float(line[2]))

    sortedHopped = {k: v for k, v in sorted(hopped.items(),reverse = True, key = lambda item:item[1][0])}
    
with open("sorted_hopped_output.txt", "w") as hopped_outfile:
    hopped_outfile.write(f"Index1-Index2\tNumber of Records\tPercent of Total Records\n")
    for stuff in sortedHopped:
        hopped_outfile.write(f"{stuff}\t{sortedHopped[stuff][0]}\t{sortedHopped[stuff][1]}\n")
    