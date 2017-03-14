"""
Average distance from a dnadist output file
"""

import os
import sys
import argparse
from roblib import stats

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='input file', required=True)
    args = parser.parse_args()

    dnadists = []
    thisdist = ""
    first = True
    with open(args.f, 'r') as f:
        for l in f:
            if first:
                # ignore the first line
                first = False
                continue
            if not l.startswith(' '):
                if thisdist:
                    dnadists.append(thisdist)
                thisdist = l.strip()
            else:
                thisdist += ' ' + l.strip()
    dnadists.append(thisdist)

    ids=[]
    pairwise = []
    for d in dnadists:
        p=d.strip().split()
        ids.append(p[0])
        pairwise.append(map(float, p[1:]))

    allpd=[]
    for i in range(len(ids)):
        for j in range(len(ids)):
            pd1 = pairwise[i][j]
            pd2 = pairwise[j][i]

            if pd1 != pd2:
                sys.stderr.write("Different diagonals in the matrix for PD1 and PD2 when i={} and j={}\n".format(i, j))
            if i == j:
                continue

            allpd.append(pd1)

    print("Average pairwise distance: {}  +/-  {}\n".format(1.0 * sum(allpd)/len(allpd), stats.stdev(allpd)))