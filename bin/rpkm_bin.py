#!/usr/bin/env python3

import sys

with open(sys.argv[1], "r") as f:
    peak = [line.strip("\n").split("\t") for line in f]

with open(sys.argv[2], "r") as f:
    bed = [line.strip("\n").split("\t") for line in f]

SIZE = int(sys.argv[3])

index = 0
n = len(peak)
num = [0] * n
for read in bed:
    mid = (int(read[1]) + int(read[2])) // 2
    while (index < n-1 and mid > int(peak[index][2])) or (index < n-1 and read[0] != peak[index][0]):
        index += 1
    num[index] += 1
    if (index < n-1) and (mid == int(peak[index+1][1])):
        num[index+1] += 1

# Output generation with formatted strings and list comprehension
SIZE_FACTOR = 10**9 / SIZE
output = [
    "{}\t{}\t{}\t{}\t{}\t{:.4f}\n".format(
        peak[i][0], peak[i][1], peak[i][2], peak[i][3], num[i],
        num[i] * SIZE_FACTOR / (int(peak[i][2]) - int(peak[i][1])) if num[i] != 0 else 0
    )
    for i in range(n)
]

with open(sys.argv[4], "w") as f:
    f.writelines(output)
