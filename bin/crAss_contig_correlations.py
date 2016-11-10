#!/usr/bin/env python

"""
Find the correlations between contigs using the crAss.pl output file output.contigs2reads.txt
"""

import os
import sys

from numpy import isnan
from scipy.stats.stats import pearsonr
import argparse

__author__ = 'Rob Edwards'


def merge_clust(c1, c2, inclust, clustermembers):
    if c2 < c1:
        [c1, c2] = [c2, c1]

    for c in clustermembers[c2]:
        clustermembers[c1].add(c)
        inclust[c] = c1

    clustermembers.pop(c2)
    return inclust, clustermembers

def notzero(x): return x > 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate pairwise pearson correlations between contigs and then cluster them')
    parser.add_argument('-d', help='data table with contigs in rows and occurence in columns', required=True)
    parser.add_argument('-s', help='minimum number of non-zero samples (default=all samples)', type=int)
    parser.add_argument('-r', help='minimum number of reads (row total) for a sample to be included (default=all rows)', type=int)
    parser.add_argument('-m', help='minimum Pearson correlation to be printed out', type=float)
    args = parser.parse_args()

    data = {}
    headers = []
    with open(args.d, 'r') as f:
        for l in f:
            p=l.strip().split("\t")
            if headers == []:
                headers = p
            else:
                tmp = list(map(int, p[1:]))
                if args.r and sum(tmp) < args.r:
                    continue
                data[p[0]] = tmp

    allcontigs = list(data.keys())
    allcontigs.sort()

    if not args.s:
        args.s = len(headers)

    # we're just going to use contigs where at least 10 samples are not zero
    nonzero = []
    for c in allcontigs:
        nz = list(filter(notzero, data[c]))
        if len(nz) >= args.s:
            nonzero.append(c)

    sys.stderr.write("Before filtering we had {} contigs, after filtering to remove samples with at least {} non zero values we have {} contigs\n".format(len(allcontigs), args.s, len(nonzero)))


    for i in range(len(nonzero)):
        cfr = nonzero[i]
        for j in range(i+1, len(nonzero)):
            cto = nonzero[j]

            dist = pearsonr(data[cfr], data[cto])[0]

            if not isnan(dist) and dist > args.m:
                print("{}\t{}\t{}".format(cfr, cto, dist))

