#!/usr/bin/env python

#script overview: takes dual_match_output.txt - made manually from summary_stats.txt and plots a 
#frequency distribution for each dual-matched index. Resulting graph to be inserted into report.md

import matplotlib.pyplot as plt
import numpy as np

indexes = []
freq = []

with open("dual_match_output.txt", "r") as dualMatched:
    dualMatched.readline() #skip header
    for line in dualMatched:
        line = line.strip()
        line = line.split("\t")
        indexes.append(line[0])
        freq.append(int(line[1]))


plt.bar(indexes, freq, color = "#008080")
plt.xticks(rotation=90)
plt.xlabel("Index")
plt.ylabel("Frequency")
plt.title("Dual Matched Index Frequency")
plt.tight_layout()
plt.savefig("dual_matched_freq.png")