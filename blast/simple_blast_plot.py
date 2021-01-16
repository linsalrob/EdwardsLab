"""
Just plot an m8 output file
"""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def run(filename, outputfile, seqname, reflen, mergey, query=False,  verbose=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if reflen == 'auto':
        minx = 100000000
        maxx = 0
        with open(filename, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                start = int(p[8])
                stop = int(p[9])
                if query:
                    if seqname and p[0] != seqname:
                        continue
                    start = int(p[6])
                    stop = int(p[7])
                elif seqname and p[1] != seqname:
                    continue
                if start > stop:
                    (start, stop) = (stop, start)
                if start < minx:
                    minx = start
                if stop > maxx:
                    maxx = stop
        maxx -= minx
        if seqname:
            sys.stderr.write(f"Auto detected the length of the sequence for {seqname} as min as {minx} max as {maxx}\n")
        else:
            sys.stderr.write(f"Auto detected the length of the sequence as min as {minx} max as {maxx}\n")
    else:
        maxx = int(reflen)
        minx = 0
        sys.stderr.write(f"Parsed min as {minx} max as {maxx}\n")

    # plot a red line for the reference
    yposn = 1
    ax.plot([minx, maxx], [yposn, yposn], color='red')
    yposn += 1


    yposns = {}
    lastyposn = yposn
    with open(filename, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            start = int(p[8])
            stop = int(p[9])
            if query:
                if seqname and p[0] != seqname:
                    continue
                start = int(p[6])
                stop = int(p[7])
            elif seqname and p[1] != seqname:
                continue
            if start > stop:
                (start, stop) = (stop, start)
            start -= minx
            stop -= minx

            if mergey:
                yidx = p[0]
                if query:
                    yidx = p[1]
                if yidx not in yposns:
                    yposns[yidx] = lastyposn
                    lastyposn += 1
                yposn = yposns[yidx]

            ax.plot([start, stop], [yposn, yposn], color='blue')
            yposn += 1

    ax.set_xlabel(f"Position in {seqname} from {minx} to {maxx}")
    ax.set_ylabel("Individual sequences")



    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.get_xaxis().tick_bottom()
    #ax.get_yaxis().tick_left()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    fig.set_facecolor('white')

    #plt.show()
    fig.savefig(outputfile)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot blast output from an m8 format file")
    parser.add_argument('-f', help='blast output filename', required=True)
    parser.add_argument('-i', help='image file to write', required=True)
    parser.add_argument('-q', help='use the query locations. The default is to use the reference locations', action='store_true')
    parser.add_argument('-n', help='name of the reference (query with -q) sequence to use. The default is to ignore this and assume you only have one sequence', default=None)
    parser.add_argument('-m', help='merge all hits with the same query (reference with -q) into a single line', action='store_true')
    parser.add_argument('-l', help='length of the reference (query with -q) sequence', default='auto')

    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    run(args.f, args.i, args.n, args.l, args.m, args.q, args.v)
