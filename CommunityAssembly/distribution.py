
"""
Given a file with the function in the first column, and then counts in
all subsequent columns, print out a two column, x-y table of the data

Add a line for a beta distribution

"""

import sys
import numpy as np
import roblib

if len(sys.argv) == 1:
    sys.exit(sys.argv[0] + " [-c number of columns as names] [-nh no header information] <output file>")

header = True
hcols = 1

for i in range(1, len(sys.argv)):
    if sys.argv[i] == "-c":
        hcols = int(sys.argv[i+1])
    if sys.argv[i] == "-nh":
        header = False

inf = sys.argv[-1]


def notzero(x): return x > 0

data = {}
total = 0
ncols = 0
with open(inf, 'r') as fin:
    if header:
        header = fin.readline()
    for l in fin:
        p = l.strip().split("\t")
        if len(p) > ncols:
            ncols = len(p)

        points = map(float, p[hcols:])
        nz = filter(notzero, points)
        psum = sum(points)
        total += psum
        data[p[0]] = [len(nz), psum]


# now calculate the mean and stdev based on the beta distribtion
xvalues = set()
betad = {}
for p in data:
    if data[p][0] not in betad:
        sys.stderr.write("alpha:" + str(data[p][0]) + " beta: " + str((ncols-data[p][0])+1)+ "\n")
        # samples = np.random.beta(data[p][0]+1, ncols-data[p][0], 1000)
        samples = np.random.beta(data[p][1]+1, total-data[p][1], 100000)
        betad[data[p][0]] = (roblib.mean(samples), roblib.stdev(samples))



seen = set()
for p in data:
    x = str(1.0 * data[p][0]/ncols)
    sys.stdout.write(p + "\t"  + str(x) + "\t" + str(1.0 * data[p][1]/total))
    if x not in seen:
        sys.stdout.write( "\t" + str(betad[data[p][0]][0]) + "\t" + str(betad[data[p][0]][1]))
    seen.add(x)
    sys.stdout.write("\n")



