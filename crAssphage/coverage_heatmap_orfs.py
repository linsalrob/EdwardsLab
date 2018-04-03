"""
Plot a heatmap of the coverage data from coverage_depth.py, except this needs the data to be transposed from that code
"""

import os
import sys

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math
matplotlib.rcParams['svg.fonttype'] = 'none'

from matplotlib.collections import PatchCollection

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot a heatmap")
    parser.add_argument('-f', help='file of tab delimited data', required=True)
    parser.add_argument('-w', help='window to use (default=1000)', default=1000, type=int)
    parser.add_argument('-n', help="column number to start at (default = 1; 0 indexed)", default=1, type=int)
    parser.add_argument('-l', help="Log normalize the data", action="store_true")
    parser.add_argument('-o', help='output file to save figure')
    parser.add_argument('-m', help='maximum value to use for any average', type=float)
    parser.add_argument('-g', help='file of THEA ORF calls', required=True)
    args = parser.parse_args()

    header = None
    data = []
    row_labels = []
    eps = 0.00000000000000001 # epsilon since we can't do log 0

    genomelength = 0
    with open(args.f, 'r') as f:
        for l in f:
            if l.startswith("#"):
                continue
            p = l.strip().split("\t")
            genomelength = len(p) - args.n
            if len(p) > args.n:
                row_labels.append(p[0])
                s = [int(x) for x in p[args.n:]]
                counter = 0
                total = 0
                thisrow = []
                for i in range(args.n, len(s)):
                    total += s[i]
                    counter += 1
                    if counter == args.w:
                        n = 1.0 * total/counter
                        if args.m and n > args.m:
                            n = args.m
                        if n:
                            if args.l:
                                thisrow.append(1.0 * math.log(n)/math.log(10))
                            else:
                                thisrow.append(n)
                        else:
                            if args.l:
                                thisrow.append(eps)
                            else:
                                thisrow.append(0)
                        counter = 0
                        total = 0
                data.append(thisrow)

    data = data[::-1]
    npd = np.array(data)

    fig, ax = plt.subplots()
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('none')

    ax.set_xlim(0, len(data[0]))
    ax.set_ylim(0, len(data))

    ax.set_xlabel("Position in genome (kbp)")
    ax.set_ylabel("Metagenome number")

    xlocs = ax.get_xticklabels()
    sys.stderr.write("There should be {} x-ticks\n".format(len(xlocs)))
    n = int((1.0 * genomelength/args.w)/len(xlocs))
    sys.stderr.write("n is: {}\n".format(n))
    sys.stderr.write("g/w is {}\n".format((1.0 * genomelength/args.w)))
    sys.stderr.write("xlocs is {}\n".format(len(xlocs)))
    xlocs = [str(int((n * i)/1000)) for i in range(len(xlocs))]
    ax.set_xticklabels(xlocs)

    heatmap = ax.pcolormesh(npd, cmap=plt.cm.Blues)
    # heatmap = ax.pcolor(allxlabelsd, npd)

    # legend
    cbar = plt.colorbar(heatmap)
    if args.l:
        cbar.set_label('log(sequence coverage)', rotation=270, labelpad=20)
    else:
        cbar.set_label('sequence coverage', rotation=270, labelpad=20)

    plt.subplots_adjust(bottom=0.3)

    pos = ax.get_position()
    # the axis go from x0 -> x1
    x0 = pos.min[0]
    x1 = pos.max[0]
    print("Axes go from {} to {}".format(x0, x1))

    # now we make a new axis for the orfs that starts at 00 and goes to 11
    # and has a transparent background
    ax2 = plt.axes([0,0,1,1], facecolor=(1,1,1,0))

    patches = []
    patch = mpatches.Rectangle([x0,0.075], x1-x0, 0.1)
    patches.append(patch)
    collection = PatchCollection(patches, color='White')
    ax2.add_collection(collection)
    patches = []

    # read in the ORF locations
    texts = []
    with open(args.g, 'r') as f:
        for l in f:
            if l.startswith('#') or l.startswith('Gene'):
                continue
            p=l.strip().split("\t")
            y = 0.125
            (start, end) = (int(p[1]), int(p[2]))
            if p[3] == '-':
                (start, end) = (int(p[2]), int(p[1]))
                y = 0.075
            start = (1.0 * start / genomelength) * (x1-x0)
            end   = (1.0 * end   / genomelength) * (x1-x0)
            width = (end - start)
            start += x0
            # sys.stderr.write("Start: {} End: {} Width: {}\n".format(start, end, width))
            # Rectangle([startx, starty], width, height
            patch = mpatches.Rectangle([start, y], width, 0.05, ec='Black', fc="DeepSkyBlue", lw=1)
            if width > 0.02:
                # we don't want to label boxes that are too small
                texts.append([start+(width/2), y + 0.025, p[0]])
            # print("{}\t{}".format(p[0], width))
            patches.append(patch)

    #p = mpatches.Rectangle([x0, 0.1], 0.5, 0.5)

    collection = PatchCollection(patches, match_original=True)
    ax2.add_collection(collection)

    for t in texts:
        ax2.text(t[0], t[1], t[2], ha='center', va='center', size='x-small')

    # plt.tight_layout()
    if args.o:
        plt.savefig(args.o)
    else:
        sys.stderr.write("Showing\n")
        plt.show()
