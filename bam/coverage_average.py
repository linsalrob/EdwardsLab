"""
Calculate the average coverage of a reference in a bam file. This just prints out a single number!
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
    parser.add_argument('-t', help='print a title that includes the name of the file (e.g. if you want to merge multiple outputs)', action='store_true')
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
        end = args.e+1

    if args.r:
        print("{}\t{}\t{}".format(args.f, args.r, (1.0 * sum(coverage[start:end]) / (end - start))))
    else:
        print("{}\t{}".format(args.f, (1.0*sum(coverage[start:end])/(end-start))))
