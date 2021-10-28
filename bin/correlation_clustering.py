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


def parse_text_file(tf, pearsoncol=2, insep="\t"):
    """
    Parse a text file and return an n-choose-2 array of the elements. The array returned has the distance from the first
    element to all other elements, and then the second element to n-1 elements (all but the first), and then the
    third element to n-2 elements (all but the first & second) and so on.
    :param tf: Text file with [a, b, pearson correlation, p-value] (e.g. the output from correlations.py)
    :type tf: str
    :param pearsoncol: the zero indexed column of the pearson correlation score
    :type numcols: int
    :param insep: Input separator. Default = tab
    :type insep: str
    :return: n-choose-2 array of the data.
    :rtype: np.array
    """

    data = {}
    ks = set()
    with open(tf, 'r') as fin:
        for l in fin:
            p=l.strip().split(insep)
            if len(p) <= pearsoncol:
                sys.stderr.write(f"ERROR: {l} does not have enough entries for {pearsoncol} to be the score. Skipped\n")
                continue
            ks.add(p[0])
            ks.add(p[1])
            if p[0] not in data:
                data[p[0]]={}
            if p[1] not in data:
                data[p[1]] = {}
            val = float(p[pearsoncol])
            if val > 1:
                val = 1
            data[p[0]][p[1]] = val
            data[p[1]][p[0]] = val

    allkeys = list(ks)
    allkeys.sort()
    # our new array will have (len(allkeys)) choose 2 elements
    nal = math.comb(len(allkeys), 2)
    nct = np.empty(nal)
    p = 0
    for i in range(len(allkeys)):
        for j in range(i + 1, len(allkeys)):
            if allkeys[j] not in  data[allkeys[i]]:
                data[allkeys[i]][allkeys[j]] = 0
            nct[p] = 1 - data[allkeys[i]][allkeys[j]]
            p += 1
    return nct


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cluster genes based on %id with cutoffs")
    parser.add_argument('-t', '--tsv', help='file with [a, b, distance] separated by tabs', required=True)
    parser.add_argument('-j', '--json', help='clusters output file name. We print them out in json format', required=True)
    parser.add_argument('-p', '--pearsoncol', help='0 indexed column in input file with the pearson score. Default = 2', type=int, default=2)
    parser.add_argument('-s', '--separator', help='Input separator. Default = tab', default="\t", type=str)
    parser.add_argument('-n', '--noclust', help='number of clusters to print (default=100)', type=int, default=100)
    args = parser.parse_args()

    matrix = parse_text_file(args.tsv, args.pearsoncol, args.separator)

    L = sch.linkage(matrix, method='average')

    print("Threshold\tNumber of clusters\tSize of largest cluster")
    with open(args.json, 'w') as out:
        out.write("[\n")
        for i in range(args.noclust+1):
            ind = sch.fcluster(L, i/args.noclust, 'distance')
            uniqs, counts = np.unique(ind, return_counts=True)
            out.write(f"{{cluster_id : {i}, largest_cluster : {max(counts)}, num_clusters: {uniqs.shape[0]}, clusters: {list(ind)}}},\n")
            print(f"{i/args.noclust}\t{uniqs.shape[0]}\t{max(counts)}")
        out.write("]\n")

