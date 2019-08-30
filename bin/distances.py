"""
Given a file with columns, generate euclidean distance between all the columns of data
"""

import os
import sys
import argparse
from itertools import combinations
from scipy.spatial import distance
import matplotlib.pyplot as plt

def read_data(f, h, c):
    """
    Read the file and store the data.
    :param f: the file to read
    :param h: whether the file contains headers
    :param c: whether the first column is a label or data
    :return:
    """

    data = []
    headers = []
    firstline = True
    start = 1
    if c:
        start = 0
    with open(f,'r') as fin:
        for l in fin:
            p=l.strip().split("\t")
            if firstline:
                if h:
                    headers = p[start:]
                for i in range(start, len(p)):
                    if not h:
                        headers.append(i)
                    data.append([])
                firstline = False
                continue
            for i in range(start, len(p)):
                data[i-start].append(float(p[i]))
    return data, headers


def pairwise(data, headers):
    """
    Calculate pairwise distances
    :param data:
    :param headers:
    :return:
    """

    cols = range(len(headers))

    for i, j in combinations(cols, 2):
        d = distance.euclidean(data[i], data[j])
        print("{}\t{}\t{}".format(headers[i], headers[j], d))


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
    parser.add_argument('-c', help='first column is data and should be included (default: the first columns is labels that are discarded)', action='store_true')
    parser.add_argument('-v', help='verbose output')
    args = parser.parse_args()

    data, headers = read_data(args.f, args.l, args.c)
    pairwise(data, headers)
    plot_pairs(data, headers)
