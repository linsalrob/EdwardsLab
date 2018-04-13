"""
Create a table of the coverage depth of reads mapping to a base in a bam file
"""

import os
import sys

import argparse
import pysam

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read a bamfile and create coverage depth (# reads that map)")
    parser.add_argument('-f', help='bamfile', required=True)
    parser.add_argument('-l', help='genome length', required=True, type=int)
    parser.add_argument('-s', help='start position (default = 0)', default=0, type=int)
    parser.add_argument('-e', help='end position (default = all)', type=int)
    parser.add_argument('-r', help="refefence name for the pileup (optional)", default=None)
    parser.add_argument('-m', help='include a column with the sum of hits after the name', action='store_true')
    parser.add_argument('-x', help='transpose the output and ommit position information. Use this if you want to combine multiple outputs', action='store_true')
    parser.add_argument('-t', help='print a title that includes the name of the file (e.g. if you want to merge multiple outputs)', action='store_true')
    parser.add_argument('-z', help='print entries with no hits. Default is to skip those', action='store_true')
    args = parser.parse_args()

    coverage = []
    for i in range(args.l+1):
        coverage.append(0)

    bam = pysam.AlignmentFile(args.f, 'rb')
    for pu in bam.pileup(reference=args.r):
        if pu.reference_pos > args.l:
            sys.stderr.write("Warning {} is larger than you genome size of {}\n".format(pu.reference_pos, args.l))
            continue
        try:
            coverage[pu.reference_pos] += pu.nsegments
        except:
            sys.stderr.write("Can't add to position {}\n".format(pu.reference_pos))

    start = 1
    end = args.l+1
    if args.s:
        start = args.s
    if args.e:
        end = args.e

    sum = 0
    for i in range(start, end):
        sum += coverage[i]

    if not args.z and 0 == sum:
        sys.exit(0)

    if args.x:
        if args.t:
            sys.stdout.write(args.f)
            if args.m:
                sys.stdout.write("\t{}".format(sum))
        else:
            # do this twice to get the tab correct
            if args.m:
                sys.stdout.write("{}".format(sum))

        for i in range(start, end):
            sys.stdout.write("\t{}".format(coverage[i]))
        sys.stdout.write("\n")
    else:
        if args.t:
            print("Position\t{}".format(args.f))
        if args.m:
            print("Sum\t{}".format(sum))
        for i in range(start, end):
            print("{}\t{}".format(i, coverage[i]))

