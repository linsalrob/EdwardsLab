"""
Given a file with the function in the first column, and then counts in
all subsequent columns, print out a two column, x-y table of the data

"""

import sys


if len(sys.argv) < 1:
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
        data[p[0]] = [psum, len(nz)]

for p in data:
    print(str(data[p][0]/total) + "\t" + str(1.0 * data[p][1]/ncols))



