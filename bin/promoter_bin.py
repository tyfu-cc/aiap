#!/usr/bin/env python3

import sys

# Read 4th column as integers from the input
with open(sys.argv[1], "r") as f:
    peak = [int(line.strip("\n").split("\t")[3]) for line in f]

# Calculate the size of each percentile bin
num = len(peak) // 100

# Compute percentile-based averages
bin = [f"{i+1}\t{sum(peak[(num*i):(num*(i+1))]) / num}\n" for i in range(99)]
bin.append(f"100\t{sum(peak[num*99:]) / num}\n")

with open(sys.argv[2], "w") as f:
    f.writelines(bin)
