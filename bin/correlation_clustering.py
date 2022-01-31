"""

Cluster correlations from correlations.py based on hierarchical linkage.

correlations.py outputs a tuple of

[object 1, object 2, pearson correlation, p value]

Here, we use 1-pearson correlation as our distance

"""
import json
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
    :return: n-choose-2 array of the data and the list of the keys in order
    :rtype: np.array, list
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
    return nct, allkeys

def generate_clusters(matrix, idlist, jsonout, noclust=100, print_singles=False):

    L = sch.linkage(matrix, method='average')

    outputdata = []
    for ele in range(noclust + 1):
        threshold = ele / noclust
        ind = sch.fcluster(L, threshold, 'distance')
        uniqs, counts = np.unique(ind, return_counts=True)
        freqs = {}
        for idx, u in enumerate(uniqs):
            freqs[u] = counts[idx]

        clusters = {}
        for idx, j in enumerate(ind):
            if freqs[j] == 1 and not print_singles:
                continue
            if j not in clusters:
                clusters[j] = []
            clusters[j].append(idlist[idx])

        # res = [i idlist[x] for x in ind]
        # res = [x for x in ind]
        # out.write(f"{{cluster_id : {i}, largest_cluster : {max(counts)}, num_clusters: {uniqs.shape[0]}, clusters: {res}}},\n")
        singles = np.isin(counts, [1]).sum()
        outputdata.append({
            "cluster_id": ele,
            "threshold": threshold,
            "largest_cluster": max(counts),
            "num_clusters": uniqs.shape[0],
            "num_singleton_clusters": singles,
            "clusters": clusters
        })
        print(f"Threshold: {threshold}\tNumber of clusters: {uniqs.shape[0]}\tLargest cluster{max(counts)}. See the output file ({jsonout}) for more details")

    with open(jsonout, 'w') as out:
        json.dump(outputdata, out)


def generate_a_cluster(matrix, idlist, jsonout, threshold=0.05, print_singles=False):
    """
    Generate a single cluster given a threshold
    """

    L = sch.linkage(matrix, method='average')

    ind = sch.fcluster(L, threshold, 'distance')
    uniqs, counts = np.unique(ind, return_counts=True)
    freqs = {}
    for idx, u in enumerate(uniqs):
        freqs[u] = counts[idx]

    clusters = {}
    for idx, j in enumerate(ind):
        if freqs[j] == 1 and not print_singles:
            continue
        ji = int(j)
        if ji not in clusters:
            clusters[ji] = []
        clusters[ji].append(idlist[idx])

    singles = int(np.isin(counts, [1]).sum())
    outputdata = {
        "threshold": threshold,
        "largest_cluster": int(max(counts)),
        "num_clusters": int(uniqs.shape[0]),
        "num_singleton_clusters": singles,
        "clusters": clusters
    }
    print(f"Threshold: {threshold}\tNumber of clusters: {uniqs.shape[0]}\tLargest cluster: {max(counts)}")

    with open(jsonout, 'w') as out:
        json.dump(outputdata, out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cluster genes based on Pearson correlation")
    parser.add_argument('-f', '--file', help='file with [a, b, distance]', required=True)
    parser.add_argument('-o', '--output', help='clusters output file name. We print them out in json format', required=True)
    parser.add_argument('-p', '--pearsoncol', help='0 indexed column in input file with the pearson score. Default = 2', type=int, default=2)
    parser.add_argument('-s', '--separator', help='Input separator. Default = tab', default="\t", type=str)
    parser.add_argument('--singles', help='print clusters with one element in them', action='store_true')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-t', '--threshold', help='clustering threshold to print a single cluster', type=float)
    group.add_argument('-n', '--numclust', help='number of clusters to print with a range of thresholds (default=100)', type=int, default=100)

    args = parser.parse_args()

    matrix, idlist = parse_text_file(args.file, args.pearsoncol, args.separator)
    if args.threshold:
        generate_a_cluster(matrix, idlist, args.output, args.threshold, args.singles)
    elif args.numclust:
        generate_clusters(matrix, idlist, args.output, args.numclust, args.singles)
