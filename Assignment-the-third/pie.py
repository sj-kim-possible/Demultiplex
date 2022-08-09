#!/usr/bin/env python

# Script overview: Plots the proportion of unknown, low quality, hopped, and dual matched records
# from demultiplexing.

import matplotlib.pyplot as plt
import numpy as np

# numbers taken from summary_stats.txt, output from demux
y = np.array([30783962, 26964891, 517612, 304980270])
mylabels = ["Unknown", "Low Quality", "Hopped", "Dual-Matched"]

colors = ['#5a7bd5','#31418b','#ffffff','#ac6ad5']

plt.pie(y, labels = mylabels, autopct='%1.1f%%', colors = colors)
plt.title(f"Percentage of Unknown, Low Quality, Hopped, and Dual-Matched")
plt.savefig("pie_results.png")
plt.show() 