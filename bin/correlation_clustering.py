"""

Cluster correlations from correlations.py based on hierarchical linkage.

correlations.py outputs a tuple of

[object 1, object 2, pearson correlation, p value]

Here, we use 1-pearson correlation as our distance

"""

import os
import sys
import math
import argparse

import numpy as np
import scipy.cluster.hierarchy as sch


def parse_text_file(tf, pearsoncol=2):
    """
    Parse a text file and return an n-choose-2 array of the elements. The array returned has the distance from the first
    element to all other elements, and then the second element to n-1 elements (all but the first), and then the
    third element to n-2 elements (all but the first & second) and so on.
    :param tf: Text file with [a, b, pearson correlation, p-value] (e.g. the output from correlations.py)
    :type tf: str
    :param pearsoncol: the zero indexed column of the pearson correlation score
    :tyoe numcols: int
    :return: n-choose-2 array of the data.
    :rtype: np.array
    """

    data = {}
    ks = set()
    with open(tf, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            if len(p) >= pearsoncol:
                sys.stderr.write(f"ERROR: {l} does not have enough entries for {pearsoncol} to be the score. Skipped\n")
                continue
            ks.add(p[0])
            ks.add(p[1])
            if p[0] not in data:
                data[p[0]]={}
            if p[1] not in data:
                data[p[1]] = {}
            data[p[0]][p[1]] = float(p[pearsoncol])
            data[p[1]][p[0]] = float(p[pearsoncol])

    allkeys = list(ks)
    allkeys.sort()
    # our new array will have (len(allkeys)) choose 2 elements
    nal = math.comb(len(allkeys), 2)
    nct = np.empty(nal)
    p = 0
    for i in range(len(allkeys)):
        for j in range(i + 1, len(allkeys)):
            nct[p] = 1 - data[allkeys[i]][allkeys[j]]
            p += 1
    return nct


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cluster genes based on %id with cutoffs")
    parser.add_argument('-t', help='file with [a, b, distance] separated by tabs', required=True)
    parser.add_argument('-j', help='clusters output file name. We print them out in json format', required=True)
    parser.add_argument('-p', help='0 indexed column in input file with the pearson score. Default = 2', type=int, default=2)
    parser.add_argument('-n', help='number of clusters to print (default=100)', type=int, default=100)
    args = parser.parse_args()

    matrix = parse_text_file(args.t, args.p)

    L = sch.linkage(matrix, method='average')

    out = open(args.j, 'w')
    for i in range(args.n+1):
        ind = sch.fcluster(L, i/args.n, 'distance')
        out.write(f"{{{i}}}: {{{ind}}},\n")
        print(f"{{{100 - i}}}\t{{{max(ind)}}}")

