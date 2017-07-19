"""
Plot a histogram of contig sizes and numbers from a directory of fasta files
"""

import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot
import numpy
from roblib import read_fasta

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', help='directory of fasta files', required=True, action='append')
    parser.add_argument('-p', help='figure file name for the graph', required=True)
    args = parser.parse_args()

    lengths = {}
    maxd = 0
    for d in args.d:
        lengths[d] = []
        for f in os.listdir(d):
            fa = read_fasta(os.path.join(d, f))
            lengths[d].extend([len(fa[x]) for x in fa])
            maxd = max(lengths[d]) if max(lengths[d]) > maxd else maxd

    bins = numpy.linspace(0, maxd, 100)
    alpha = 1.0 / len(args.d)

    for d in args.d:
        pyplot.hist(lengths[d], bins, alpha=alpha, label=d)
    pyplot.legend(loc='upper right')
    pyplot.savefig(args.p)
