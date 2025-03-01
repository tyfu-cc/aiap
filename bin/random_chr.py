#!/usr/bin/env python3

import sys

batch_size = 1e4
output = []

with open(sys.argv[1], "r") as infile, open(sys.argv[2], "w") as outfile:
    for line in infile:
        line = line.strip().split("\t")
        chr_name, chr_size = line[0], int(line[1])

        # Calculate the range in chunks of 500
        for i in range(chr_size // 500):
            output.append("{0}\t{1}\t{2}\trandom_non_peak_{3}\n".format(chr_name, i * 500, (i + 1) * 500, i + 1))
            # outfile.write(chr_name + "\t" + str(i * 500) + "\t" + str((i + 1) * 500) + "\t" + "random_non_peak_" + str(i + 1) + "\n")
            if len(output) >= batch_size:
                outfile.writelines(output)
                output = []
    if output:
        outfile.writelines(output)
