"""
Plot coverage of the three primer regions from a given bam file
"""

import os
import sys
import argparse

import pysam
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy


def get_coverage(bamfile, genomelen, refg):
    """
    Get the coverage of the bamfile from 0 to genomelen. Returns an array with each element the coverage.

    :param bamfile: the bam file to read
    :param genomelen: the length of the genome
    :param refg: the reference genome
    :return: an array of coverage
    """

    coverage = []
    for i in range(genomelen):
        coverage.append(0)

    bam = pysam.AlignmentFile(bamfile, 'rb')
    for pu in bam.pileup(reference=refg):
        try:
            coverage[pu.reference_pos] += pu.nsegments
        except:
            sys.stderr.write("Can't add to position {}\n".format(pu.reference_pos))

    return coverage


def plotcoverage(x, y, plotax, ymax, linecolor, xtitle, ytitle):
    """
    Plot the coverage. x and y are two arrays of the same length
    :param x: x-coordinates
    :param y: y-coordinates
    :param plotax: the axes to plot upon
    :param ymax: the max value for the y axis
    :param linecolor: the color of the line
    :param xtitle: the title for the x axis
    :param ytitle: the title for the y axis
    :return:
    """

    if len(x) != len(y):
        sys.stderr.write("The x-array length ({}) and y-array length ({}) must be the same\n".format(len(x), len(y)))
        sys.exit(-1)

    plotax.set_xlabel(xtitle)
    plotax.set_ylabel(ytitle)
    plotax.set_ylim(0,ymax)
    plotax.plot(x, y, color=linecolor, marker='None', linestyle='-')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', help='bam file', default='/home/redwards/Desktop/GlobalCrassphage/gretel_reruns/SRR4408002/SRR4408002.bam')
    parser.add_argument('-o', help='plot on one figure panel (default is 3 panes)', action='store_true')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    fig = plt.figure()

    # our three primer regions
    locations = {
            'PrimerA': (25633, 26964),
            'PrimerB': (33708, 35062),
            'PrimerC': (43819, 45057)
    }

    # retrieve all coverage
    coverage = get_coverage(args.b, 97092, 'JQ995537')
    ymaxes=[]
    for i, primer in enumerate(['PrimerA', 'PrimerB', 'PrimerC']):
        ymaxes.append(max(coverage[locations[primer][0]:locations[primer][1]+1]))
    ymax = max(ymaxes) # used to scale the y axis
    ymax += int(ymax/10)

    if args.o:
        ax = fig.add_subplot(111)
        linecolors=['blue', 'red', 'green']
        legends = []
        for i, primer in enumerate(['PrimerA', 'PrimerB', 'PrimerC']):
            x = range(locations[primer][1]-locations[primer][0]+1)
            y = coverage[locations[primer][0]:locations[primer][1]+1]
            plotcoverage(x, y, ax, ymax, linecolors[i], "Coverage", "Average Read Depth")
            legends.append(ax)
        ax.legend(legends, ['PrimerA', 'PrimerB', 'PrimerC'])
    else:
        for i, primer in enumerate(['PrimerA', 'PrimerB', 'PrimerC']):
            ax = fig.add_subplot(3, 1, i + 1)
            x = range(locations[primer][0], locations[primer][1] + 1)
            y = coverage[locations[primer][0]:locations[primer][1] + 1]
            plotcoverage(x, y, ax, ymax, "blue", "Coverage of {}".format(primer), "Average Read Depth")

    fig.set_facecolor('white')

    plt.tight_layout()
    plt.show()
    fig.savefig('/home/redwards/Desktop/GlobalCrassphage/gretel_reruns/SRR4408002/SRR4408002.png')