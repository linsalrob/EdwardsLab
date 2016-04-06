"""

Cluster genes based on %id with cutoffs

"""

import os
import sys

import argparse

import scipy
import scipy.cluster.hierarchy as sch

def parse_text_file(tf):
    """
    Parse a text file and return an n-choose-2 array of the elements. The array returned has the distance from the first
    element to all other elements, and then the second element to n-1 elements (all but the first), and then the
    third element to n-2 elements (all but the first & second) and so on.
    :param tf: Text file with [a, b, distance]
    :type tf: str
    :return: n-choose-2 array of the data.
    :rtype: array
    """

    data = {}
    ks = set()
    with open(tf, 'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            ks.add(p[0])
            ks.add(p[1])
            if p[0] not in data:
                data[p[0]]={}
            if p[1] not in data:
                data[p[1]] = {}
            data[p[0]][p[1]] = float(p[2])/100
            data[p[1]][p[0]] = float(p[2])/100

    allkeys = list(ks)
    allkeys.sort()
    nct = []
    for i in range(len(allkeys)):
        for j in range(i+1, len(allkeys)):
            nct.append(data[allkeys[i]][allkeys[j]])
    return nct

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cluster genes based on %id with cutoffs")
    parser.add_argument('-t', help='file with [a, b, distance] separated by tabs', required=True)
    args = parser.parse_args()

    matrix = parse_text_file(args.t)
    L = sch.linkage(matrix, method='complete')

    for i in range(101):
        ind = sch.fcluster(L, i/100.0, 'distance')
        print("{}\t{}".format(100-i, max(ind)))

