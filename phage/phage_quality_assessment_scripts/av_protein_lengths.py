"""
Calculate the average protein lengths from the blast output files.
"""

import os
import sys
import argparse
from roblib import bcolors, median

def av_protein_lengths(sample, blastfile, fractionout, summaryout, searchtype):
    """
    Calculate the average length of the best hit of all the proteins
    """
    q = {}
    av = []
    sys.stderr.write(f"{bcolors.GREEN}Average protein lengths for {sample} and {searchtype}{bcolors.ENDC}\n")

    with open(blastfile, 'r') as f:
        with open(fractionout, 'w') as out:
            for l in f:
                p = l.strip().split("\t")
                if p[0] in q:
                    continue
                q[p[0]] = int(p[12])/int(p[13])
                av.append(q[p[0]])
                out.write(f"{p[0]}\t{q[p[0]]}\n")
    with open(summaryout, 'w') as out:
        out.write(f"{sample}\tAverage {searchtype} protein lengths\t")
        out.write("[num orfs, median proportional length, average proportional length]\t")
        out.write(f"{len(av)}\t{median(av)}\t{sum(av)/len(av)}\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-s', help='sample name used in output', required=True)
    parser.add_argument('-b', help='blast m8 file', required=True)
    parser.add_argument('-f', help='fractions output file', required=True)
    parser.add_argument('-o', help='summary output file', required=True)
    parser.add_argument('-t', help='search type  (e.g. phage, bacteria) (used in output)', required=True)
    args = parser.parse_args()

    av_protein_lengths(args.s, args.b, args.f, args.o, args.t)