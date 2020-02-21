"""
Calculate the kurtosis of a bam file. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4424905/ for the definition
"""

import os
import sys

import argparse
import pysam

from roblib import mean, stdev

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read a bamfile calculate average coverage and kurtosis")
    parser.add_argument('-f', help='bamfile', required=True)
    parser.add_argument('-l', help='genome length', required=True, type=int)
    parser.add_argument('-s', help='start position (default = 0)', default=0, type=int)
    parser.add_argument('-e', help='end position (default = all)', type=int)
    parser.add_argument('-r', help="refefence name for the pileup (optional)", default=None)
    parser.add_argument('-c', help='print column names for the output', action='store_true')
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

    # here we trim the coverage array to make the math easier!
    start = 1
    end = args.l+1
    if args.e:
        end = args.e+1
        coverage = coverage[0:end]
    if args.s:
        start = args.s
        coverage = coverage[start:]

    # calculate the average over coverage
    av = mean(coverage)
    st = stdev(coverage)
    
    k = 0
    for i,j in enumerate(coverage):
        k += (j-av)**4
    
    k = k / (len(coverage) * (st**4))

    k -= 3


    if args.r:
        if args.c:
            print("Filename\tReference\tAverage\tStDev\tKurtosis")
        print(f"{args.f}\t{args.r}\t{av}\t{st}\t{k}")
    else:
        if args.c:
            print("Filename\tAverage\tStDev\tKurtosis")
        print(f"{args.f}\t{av}\t{st}\t{k}")
