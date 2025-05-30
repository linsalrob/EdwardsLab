#!/usr/bin/env python3

"""
Generate correlations between all the columns or rows of data. We will read a matrix file
and transpose it if requested, to generate all pairwise Pearson's correlations.

Note see correlation_clustering.py to cluster these using hierarchical linkage
"""

import os
import sys
import argparse
from itertools import combinations
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

def read_data(f, h, c, verbose=False):
    """
    Read the file and store the data.
    :param f: the file to read
    :param h: whether the file contains headers
    :param c: whether the first column is a label or data
    :return:
    """

    if verbose:
        sys.stderr.write(f"Reading data from {f}. Column headers: {h} Row names: {~c}\n")

    data = []
    headers = []
    rownames = []
    firstline = True
    start = 1
    if c:
        start = 0
    n = 0
    with open(f,'r') as fin:
        for l in fin:
            p=l.rstrip().split("\t")
            if firstline:
                if h:
                    headers = p[start:]
                for i in range(start, len(p)):
                    if not h:
                        headers.append(i)
                    data.append([])
                firstline = False
                continue
            if start == 1:
                rownames.append(p[0])
            else:
                n+=1
                rownames.append(f"row{n}")
            for i in range(start, len(p)):
                data[i-start].append(float(p[i]))
    return data, headers, rownames


def transpose(data):
    """
    Transpose the data. In case the data are in rows, not columns
    """
    return [*zip(*data)]

def pairwise(data, headers, contigname=None):
    """
    Calculate pairwise distances
    :param data: the raw data matrix
    :param headers: the array of contig names
    :return:
    """

    cols = range(len(headers))

    for i, j in combinations(cols, 2):
        if contigname and (headers[i] != contigname and headers[j] != contigname):
            continue
        pearson, p = pearsonr(data[i], data[j])
        print("{}\t{}\t{}\t{}".format(headers[i], headers[j], pearson, p))


def plot_pairs(data, headers):
    cols = range(len(headers))
    f, axarr = plt.subplots(2, 2)

    pltx = 0
    plty = 0
    for i, j in combinations(cols, 2):
        axarr[pltx, plty].plot(data[i], data[j], 'ro')
        axarr[pltx, plty].set_title('{} versus {}'.format(headers[i], headers[j]))
        pltx += 1
        if pltx == 2:
            pltx = 0
            plty = 1

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate all pairwise correlations between the data")
    parser.add_argument('-f', help='file of data with data in columns', required=True)
    parser.add_argument('-l', help='first line is a header line (will be used in output)', action='store_true')
    parser.add_argument('-c', help='first column is data and should be included (default: the first columns is labels that are saved if transposed)', action='store_true')
    parser.add_argument('-r', help='data is in rows (ie. transpose)', action='store_true')
    parser.add_argument('-n', help='filter for a particular contig that must be included in the pairwise calculations')
    parser.add_argument('-p', help='plot the pairwise data', action='store_true')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    data, headers, rownames = read_data(args.f, args.l, args.c)
    if args.r:
        data = transpose(data)
        headers = rownames

    pairwise(data, headers, args.n)
    if args.p:
        plot_pairs(data, headers)
